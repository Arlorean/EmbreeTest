// ======================================================================== //
// Copyright 2009-2017 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#include "<embree3/rtcore.h>"
#include "<embree3/rtcore_ray.h>"

#include <ppl.h>
using namespace concurrency;

#include <malloc.h>
#define alignedMalloc(s) _aligned_malloc(s, 64)
#define alignedFree(p) _aligned_free(p)
#define __aligned(n) __declspec(align(n))

#include "linalg.h"
using namespace linalg;
using namespace linalg::aliases;

#include <iostream>
#include <fstream>
#include <chrono>

namespace linalg
{
	// https://msdn.microsoft.com/en-us/library/windows/desktop/bb281710(v=vs.85).aspx
	template <class T>
	mat<T, 4, 4> lookat_matrix(const vec<T, 3> &eye, const vec<T, 3> &target, vec<T, 3> &up)
	{
		auto zaxis = normalize(target - eye);
		auto xaxis = normalize(cross(up, zaxis));
		auto yaxis = cross(zaxis, xaxis);
		yaxis = normalize(yaxis);
		auto orientation = mat<T, 4, 4>(
			vec<T, 4>(xaxis.x, yaxis.x, zaxis.x, 0),
			vec<T, 4>(xaxis.y, yaxis.y, zaxis.y, 0),
			vec<T, 4>(xaxis.z, yaxis.z, zaxis.z, 0),
			vec<T, 4>(0, 0, 0, 1));
		auto translation = translation_matrix(-eye);
		return mul(orientation, translation);
	}

	template <class T>
	vec<T, 3> mul(const mat<T, 4, 4> &m, const vec<T, 3> &v)
	{
		auto v4 = vec<T, 4>(v.x, v.y, v.z, 1);
		return mul(m, v4).xyz();
	}
} // namespace linalg

namespace embree
{
	struct RTC_ALIGN(32) RTCValid8
	{
		int data[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
	};

	struct Color
	{
		char r, g, b;
	};
	struct Vertex4f
	{
		float x, y, z, r;
	};
	struct Triangle
	{
		int v0, v1, v2;
	};

	inline float clamp(float val, float min, float max) { return (val < min ? min : (val > max ? max : val)); }
	inline Color color(float3 v)
	{
		Color c;
		c.r = (char)(255.0f * clamp(v.x, 0, 1));
		c.g = (char)(255.0f * clamp(v.y, 0.0f, 1.0f));
		c.b = (char)(255.0f * clamp(v.z, 0.0f, 1.0f));
		return c;
	}

	/* scene data */
	RTCDevice g_device = nullptr;
	RTCScene g_scene = nullptr;
	float3 *face_colors = nullptr;
	float3 *vertex_colors = nullptr;

	/* adds a cube to the scene */
	unsigned int addCube(RTCScene scene_i)
	{
		/* create a triangulated cube with 12 triangles and 8 vertices */
		RTCGeometry mesh = rtcNewGeometry(g_device, RTC_GEOMETRY_TYPE_TRIANGLE);

		/* create face and vertex color arrays */
		face_colors = (float3 *)alignedMalloc(12 * sizeof(float3), 16);
		vertex_colors = (float3 *)alignedMalloc(8 * sizeof(float3), 16);

		/* set vertices and vertex colors */
		Vertex4f *vertices = (Vertex4f *)rtcSetNewGeometryBuffer(mesh, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, sizeof(Vertex4f), 8);
		vertex_colors[0] = float3(0, 0, 0);
		vertices[0].x = -1;
		vertices[0].y = -1;
		vertices[0].z = -1;
		vertex_colors[1] = float3(0, 0, 1);
		vertices[1].x = -1;
		vertices[1].y = -1;
		vertices[1].z = +1;
		vertex_colors[2] = float3(0, 1, 0);
		vertices[2].x = -1;
		vertices[2].y = +1;
		vertices[2].z = -1;
		vertex_colors[3] = float3(0, 1, 1);
		vertices[3].x = -1;
		vertices[3].y = +1;
		vertices[3].z = +1;
		vertex_colors[4] = float3(1, 0, 0);
		vertices[4].x = +1;
		vertices[4].y = -1;
		vertices[4].z = -1;
		vertex_colors[5] = float3(1, 0, 1);
		vertices[5].x = +1;
		vertices[5].y = -1;
		vertices[5].z = +1;
		vertex_colors[6] = float3(1, 1, 0);
		vertices[6].x = +1;
		vertices[6].y = +1;
		vertices[6].z = -1;
		vertex_colors[7] = float3(1, 1, 1);
		vertices[7].x = +1;
		vertices[7].y = +1;
		vertices[7].z = +1;

		/* set triangles and face colors */
		int tri = 0;
		Triangle *triangles = (Triangle *)rtcSetNewGeometryBuffer(mesh, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, sizeof(Triangle), 12);

		// left side
		face_colors[tri] = float3(1, 0, 0);
		triangles[tri].v0 = 0;
		triangles[tri].v1 = 1;
		triangles[tri].v2 = 2;
		tri++;
		face_colors[tri] = float3(1, 0, 0);
		triangles[tri].v0 = 1;
		triangles[tri].v1 = 3;
		triangles[tri].v2 = 2;
		tri++;

		// right side
		face_colors[tri] = float3(0, 1, 0);
		triangles[tri].v0 = 4;
		triangles[tri].v1 = 6;
		triangles[tri].v2 = 5;
		tri++;
		face_colors[tri] = float3(0, 1, 0);
		triangles[tri].v0 = 5;
		triangles[tri].v1 = 6;
		triangles[tri].v2 = 7;
		tri++;

		// bottom side
		face_colors[tri] = float3(0.5f);
		triangles[tri].v0 = 0;
		triangles[tri].v1 = 4;
		triangles[tri].v2 = 1;
		tri++;
		face_colors[tri] = float3(0.5f);
		triangles[tri].v0 = 1;
		triangles[tri].v1 = 4;
		triangles[tri].v2 = 5;
		tri++;

		// top side
		face_colors[tri] = float3(1.0f);
		triangles[tri].v0 = 2;
		triangles[tri].v1 = 3;
		triangles[tri].v2 = 6;
		tri++;
		face_colors[tri] = float3(0.5f, 0.0f, 0.0f);
		triangles[tri].v0 = 3;
		triangles[tri].v1 = 7;
		triangles[tri].v2 = 6;
		tri++;

		// front side
		face_colors[tri] = float3(0, 0, 1);
		triangles[tri].v0 = 0;
		triangles[tri].v1 = 2;
		triangles[tri].v2 = 4;
		tri++;
		face_colors[tri] = float3(0, 0, 1);
		triangles[tri].v0 = 2;
		triangles[tri].v1 = 6;
		triangles[tri].v2 = 4;
		tri++;

		// back side
		face_colors[tri] = float3(1, 1, 0);
		triangles[tri].v0 = 1;
		triangles[tri].v1 = 5;
		triangles[tri].v2 = 3;
		tri++;
		face_colors[tri] = float3(1, 1, 0);
		triangles[tri].v0 = 3;
		triangles[tri].v1 = 5;
		triangles[tri].v2 = 7;
		tri++;

		rtcSetGeometryVertexAttributeCount(mesh, 1);
		rtcSetSharedGeometryBuffer(mesh, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, RTC_FORMAT_FLOAT3, vertex_colors, 0, sizeof(float3), 8);

		rtcCommitGeometry(mesh);
		unsigned int geomID = rtcAttachGeometry(scene_i, mesh);
		rtcReleaseGeometry(mesh);

		return geomID;
	}

	/* adds a ground plane to the scene */
	unsigned int addGroundPlane(RTCScene scene_i)
	{
		/* create a triangulated plane with 2 triangles and 4 vertices */
		RTCGeometry mesh = rtcNewGeometry(g_device, RTC_GEOMETRY_TYPE_TRIANGLE);

		/* set vertices */
		Vertex4f *vertices = (Vertex4f *)rtcSetNewGeometryBuffer(mesh, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, sizeof(Vertex4f), 4);
		vertices[0].x = -10;
		vertices[0].y = -2;
		vertices[0].z = -10;
		vertices[1].x = -10;
		vertices[1].y = -2;
		vertices[1].z = +10;
		vertices[2].x = +10;
		vertices[2].y = -2;
		vertices[2].z = -10;
		vertices[3].x = +10;
		vertices[3].y = -2;
		vertices[3].z = +10;

		/* set triangles */
		Triangle *triangles = (Triangle *)rtcSetNewGeometryBuffer(mesh, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, sizeof(Triangle), 2);
		triangles[0].v0 = 0;
		triangles[0].v1 = 1;
		triangles[0].v2 = 2;
		triangles[1].v0 = 1;
		triangles[1].v1 = 3;
		triangles[1].v2 = 2;

		rtcCommitGeometry(mesh);
		unsigned int geomID = rtcAttachGeometry(scene_i, mesh);
		rtcReleaseGeometry(mesh);
		return geomID;
	}

	static const float3 lightDir = normalize(float3(-1, -1, -1));

	float3 castRay(float3 rayOrg, float3 rayDir)
	{
		/* initialize ray */
		RTCRayHit ray;
		ray.ray.flags = 0;

		memcpy(&ray.ray.org_x, &rayOrg.x, sizeof(float));
		memcpy(&ray.ray.org_y, &rayOrg.y, sizeof(float));
		memcpy(&ray.ray.org_z, &rayOrg.z, sizeof(float));
		memcpy(&ray.ray.dir_x, &rayDir.x, sizeof(float));
		memcpy(&ray.ray.dir_y, &rayDir.y, sizeof(float));
		memcpy(&ray.ray.dir_z, &rayDir.z, sizeof(float));
		ray.ray.tnear = 0.001f;
		ray.ray.tfar = INFINITY;
		ray.hit.geomID = RTC_INVALID_GEOMETRY_ID;
		ray.hit.primID = RTC_INVALID_GEOMETRY_ID;
		ray.ray.mask = -1;
		ray.ray.time = 0;

		/* intersect ray with scene */
		{
			RTCIntersectContext context;
			rtcInitIntersectContext(&context);
			rtcIntersect1(g_scene, &context, &ray);
			ray.hit.Ng_x = -ray.hit.Ng_x;
			ray.hit.Ng_y = -ray.hit.Ng_y;
			ray.hit.Ng_z = -ray.hit.Ng_z;
		}

		/* shade pixels */
		float3 color = float3(0.0f);
		if (ray.hit.geomID != RTC_INVALID_GEOMETRY_ID)
		{
			float3 diffuse = face_colors[ray.hit.primID];
			color = color + diffuse * 0.5f;

			/* initialize shadow ray */
			auto shadowOrg = (rayOrg + rayDir * ray.ray.tfar);
			auto shadowDir = (lightDir * -1.0f);
			RTCRay shadow;
			shadow.flags = 0;
			memcpy(&shadow.org_x, &shadowOrg.x, sizeof(float));
			memcpy(&shadow.org_y, &shadowOrg.y, sizeof(float));
			memcpy(&shadow.org_z, &shadowOrg.z, sizeof(float));
			memcpy(&shadow.dir_x, &shadowDir.x, sizeof(float));
			memcpy(&shadow.dir_y, &shadowDir.y, sizeof(float));
			memcpy(&shadow.dir_z, &shadowDir.z, sizeof(float));
			shadow.tnear = 0.001f;
			shadow.tfar = INFINITY;
			shadow.mask = -1;
			shadow.time = 0;

			/* trace shadow ray */
			{
				RTCIntersectContext context;
				rtcInitIntersectContext(&context);
				rtcOccluded1(g_scene, &context, &shadow);
				// EMBREE_FIXME: shadow is occluded when shadow.ray.tfar < 0.0f
			}

			/* add light contribution */
			if (shadow.tfar >= 0.0f)
			{
				color = color + diffuse * clamp(-dot(lightDir, normalize(float3(ray.hit.Ng_x, ray.hit.Ng_y, ray.hit.Ng_z))), 0.0f, 1.0f);
			}
		}
		return color;
	}

	void castRay8(const float3 rayOrg, const float3 *rayDir, float3 *colors)
	{
		/* initialize ray */
		RTCRayHit8 ray;
		for (auto i = 0; i < 8; ++i)
		{
			ray.ray.org_x[i] = rayOrg.x;
			ray.ray.org_y[i] = rayOrg.y;
			ray.ray.org_z[i] = rayOrg.z;
			ray.ray.dir_x[i] = rayDir[i].x;
			ray.ray.dir_y[i] = rayDir[i].y;
			ray.ray.dir_z[i] = rayDir[i].z;

			ray.ray.tnear[i] = 0.0001f;
			ray.ray.tfar[i] = INFINITY;
			ray.hit.geomID[i] = RTC_INVALID_GEOMETRY_ID;
			ray.hit.primID[i] = RTC_INVALID_GEOMETRY_ID;
			ray.ray.mask[i] = -1;
			ray.ray.time[i] = 0;
		}

		/* intersect ray with scene */
		RTCValid8 valid;
		{
			RTCIntersectContext context;
			rtcInitIntersectContext(&context);
			rtcIntersect8(valid.data, g_scene, &context, &ray);
		}

		/* initialize shadow ray */
		RTCRay8 shadow;
		for (auto i = 0; i < 8; ++i)
		{
			valid.data[i] = (ray.hit.geomID[i] != RTC_INVALID_GEOMETRY_ID) ? -1 : 0;
			if (valid.data[i])
			{
				auto shadowOrg = (rayOrg + rayDir[i] * ray.ray.tfar[i]);
				auto shadowDir = (lightDir * -1.0f);
				shadow.org_x[i] = shadowOrg.x;
				shadow.org_y[i] = shadowOrg.y;
				shadow.org_z[i] = shadowOrg.z;
				shadow.dir_x[i] = shadowDir.x;
				shadow.dir_y[i] = shadowDir.y;
				shadow.dir_z[i] = shadowDir.z;
				shadow.tnear[i] = 0.001f;
				shadow.tfar[i] = INFINITY;
				shadow.mask[i] = -1;
				shadow.time[i] = 0;
			}
		}

		/* trace shadow ray */
		{
			RTCIntersectContext context;
			rtcInitIntersectContext(&context);
			rtcOccluded8(valid.data, g_scene, &context, &shadow);
		}

		/* shade pixels */
		for (auto i = 0; i < 8; ++i)
		{
			float3 color = float3(0.0f);

			if (valid.data[i])
			{
				float3 diffuse = face_colors[ray.hit.primID[i]];
				color = color + diffuse * 0.5f;

				/* add light contribution */
				if (shadow.tfar[i] >= 0)
				{
					color = color + diffuse * clamp(-dot(lightDir, normalize(float3(ray.hit.Ng_x[i], ray.hit.Ng_y[i], ray.hit.Ng_z[i]))), 0.0f, 1.0f);
				}
			}

			colors[i] = color;
		}
	}
} // namespace embree

using namespace embree;
using namespace std;

void ViewImageFile(const wstring &fileName);

void main(int argc, char *argv[])
{
	// Image
	auto width = 512, height = 512;
	auto pixels = new Color[width * height];

	// Embree device and scene
	g_device = rtcNewDevice(nullptr);
	g_scene = rtcNewScene(g_device);
	rtcSetSceneFlags(g_scene, RTC_SCENE_FLAG_COMPACT);
	rtcSetSceneBuildQuality(g_scene, RTC_BUILD_QUALITY_MEDIUM);
	addCube(g_scene);
	addGroundPlane(g_scene);
	rtcCommitScene(g_scene);

	// Camera
	auto from = float3(1.5, 1.5, -1.5);
	auto to = float3();
	auto up = float3(0, 1, 0);
	auto viewMatrix = lookat_matrix(from, to, up);
	auto viewMatrixInverse = inverse(viewMatrix);

	// Projection
	auto fovy = (float)(90.0 / 180 * 3.141592653589); // in radians
	auto z = 1 / tan(fovy / 2.0f);
	auto aspect = (float)width / height;
	auto dx = 1.0f / width, dy = 1.0f / height;

	// Render scene using single rays
	{
		auto startTime = chrono::high_resolution_clock::now();
		parallel_for(size_t(0), size_t(width * height), [&](size_t i) {
			auto x = (i % width);
			auto y = (i / width);
			auto u = (2 * (x * dx) - 1);
			auto v = (2 * (y * dy) - 1);
			auto d = mul(viewMatrixInverse, float3(u * aspect, -v, z)) - from;
			auto c = castRay(from, normalize(d));
			pixels[i] = color(c);
		});
		auto endTime = chrono::high_resolution_clock::now();
		cout << "Render time (1) = " << chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() << endl;
	}

	// Render scene using 8 ray packets
	{
		auto startTime = chrono::high_resolution_clock::now();
		parallel_for(size_t(0), size_t(width * height / 8), [&](size_t _i) {
			float3 rays[8];
			for (auto r = 0; r < 8; ++r)
			{
				auto i = (_i * 8) + r;
				auto x = (i % width);
				auto y = (i / width);
				auto u = (2 * (x * dx) - 1);
				auto v = (2 * (y * dy) - 1);
				auto d = mul(viewMatrixInverse, float3(u * aspect, -v, z)) - from;
				rays[r] = normalize(d);
			}
			float3 colors[8];
			castRay8(from, rays, colors);
			for (auto r = 0; r < 8; ++r)
			{
				auto i = (_i * 8) + r;
				pixels[i] = color(colors[r]);
			}
		});
		auto endTime = chrono::high_resolution_clock::now();
		cout << "Render time (8) = " << chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() << endl;
	}

	// Free device and scene
	rtcReleaseScene(g_scene);
	g_scene = nullptr;
	rtcReleaseDevice(g_device);
	g_device = nullptr;
	alignedFree(face_colors);
	face_colors = nullptr;
	alignedFree(vertex_colors);
	vertex_colors = nullptr;

	// Write image file
	auto file = ofstream(L"image.ppm");
	file << "P3" << endl;
	file << width << " " << height << endl;
	file << 255 << endl;
	for (int i = 0; i < width * height; i++)
	{
		auto c = pixels[i];
		file << (int)pixels[i].r << " " << (int)pixels[i].g << " " << (int)pixels[i].b << " \n";
	}
	file.close();

	// Free image
	delete[] pixels;

	// Open image with default viewer
	ViewImageFile(L"image.ppm");
}

#include <windows.h>

void ViewImageFile(const wstring &fileName)
{
	ShellExecute(0, 0, fileName.c_str(), 0, 0, SW_SHOW);
}
