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

#include <embree2/rtcore.h>
#include <embree2/rtcore_ray.h>

#include <ppl.h>
using namespace concurrency;

#include <malloc.h>  
#define alignedMalloc(s) _aligned_malloc(s,64)
#define alignedFree(p) _aligned_free(p)
#define __aligned(n) __declspec(align(n))

#include "linalg.h"
using namespace linalg;
using namespace linalg::aliases;

#include <iostream>
#include <fstream>
#include <chrono>

namespace linalg {
	// https://msdn.microsoft.com/en-us/library/windows/desktop/bb281710(v=vs.85).aspx
	template<class T> mat<T, 4, 4> lookat_matrix(const vec<T,3>& eye, const vec<T, 3>& target, vec<T, 3>& up) {
		auto zaxis = normalize(target-eye);
		auto xaxis = normalize(cross(up, zaxis));
		auto yaxis = cross(zaxis, xaxis);
		yaxis = normalize(yaxis);
		auto orientation = mat<T, 4, 4>(
			vec<T,4>(xaxis.x, yaxis.x, zaxis.x, 0),
			vec<T,4>(xaxis.y, yaxis.y, zaxis.y, 0),
			vec<T,4>(xaxis.z, yaxis.z, zaxis.z, 0),
			vec<T,4>(0,0,0,1)
		);
		auto translation = translation_matrix(-eye);
		return mul(orientation, translation);
	}

	template<class T> vec<T, 3> mul(const mat<T, 4, 4>& m, const vec<T, 3>& v) {
		auto v4 = vec<T,4>(v.x, v.y, v.z, 1);
		return mul(m, v4).xyz();
	}
}

namespace embree {

	struct Color { char r, g, b; };
	struct Vertex { float x, y, z, r; }; // FIXME: rename to Vertex4f
	struct Triangle { int v0, v1, v2; };

	inline float clamp(float val, float min, float max) { return (val < min ? min : (val > max ? max : val)); }
	inline Color color(float3 v) {
		Color c;
		c.r = (char)(255.0f * clamp(v.x, 0, 1));
		c.g = (char)(255.0f * clamp(v.y, 0.0f, 1.0f));
		c.b = (char)(255.0f * clamp(v.z, 0.0f, 1.0f));
		return c;
	}

	/* scene data */
	RTCDevice g_device = nullptr;
	RTCScene g_scene = nullptr;
	float3* face_colors = nullptr;
	float3* vertex_colors = nullptr;

	/* adds a cube to the scene */
	unsigned int addCube(RTCScene scene_i)
	{
		/* create a triangulated cube with 12 triangles and 8 vertices */
		unsigned int mesh = rtcNewTriangleMesh(scene_i, RTC_GEOMETRY_STATIC, 12, 8);

		/* create face and vertex color arrays */
		face_colors = (float3*)alignedMalloc(12 * sizeof(float3));
		vertex_colors = (float3*)alignedMalloc(8 * sizeof(float3));

		/* set vertices and vertex colors */
		Vertex* vertices = (Vertex*)rtcMapBuffer(scene_i, mesh, RTC_VERTEX_BUFFER);
		vertex_colors[0] = float3(0, 0, 0); vertices[0].x = -1; vertices[0].y = -1; vertices[0].z = -1;
		vertex_colors[1] = float3(0, 0, 1); vertices[1].x = -1; vertices[1].y = -1; vertices[1].z = +1;
		vertex_colors[2] = float3(0, 1, 0); vertices[2].x = -1; vertices[2].y = +1; vertices[2].z = -1;
		vertex_colors[3] = float3(0, 1, 1); vertices[3].x = -1; vertices[3].y = +1; vertices[3].z = +1;
		vertex_colors[4] = float3(1, 0, 0); vertices[4].x = +1; vertices[4].y = -1; vertices[4].z = -1;
		vertex_colors[5] = float3(1, 0, 1); vertices[5].x = +1; vertices[5].y = -1; vertices[5].z = +1;
		vertex_colors[6] = float3(1, 1, 0); vertices[6].x = +1; vertices[6].y = +1; vertices[6].z = -1;
		vertex_colors[7] = float3(1, 1, 1); vertices[7].x = +1; vertices[7].y = +1; vertices[7].z = +1;
		rtcUnmapBuffer(scene_i, mesh, RTC_VERTEX_BUFFER);

		/* set triangles and face colors */
		int tri = 0;
		Triangle* triangles = (Triangle*)rtcMapBuffer(scene_i, mesh, RTC_INDEX_BUFFER);

		// left side
		face_colors[tri] = float3(1, 0, 0); triangles[tri].v0 = 0; triangles[tri].v1 = 2; triangles[tri].v2 = 1; tri++;
		face_colors[tri] = float3(1, 0, 0); triangles[tri].v0 = 1; triangles[tri].v1 = 2; triangles[tri].v2 = 3; tri++;

		// right side
		face_colors[tri] = float3(0, 1, 0); triangles[tri].v0 = 4; triangles[tri].v1 = 5; triangles[tri].v2 = 6; tri++;
		face_colors[tri] = float3(0, 1, 0); triangles[tri].v0 = 5; triangles[tri].v1 = 7; triangles[tri].v2 = 6; tri++;

		// bottom side
		face_colors[tri] = float3(0.5f);  triangles[tri].v0 = 0; triangles[tri].v1 = 1; triangles[tri].v2 = 4; tri++;
		face_colors[tri] = float3(0.5f);  triangles[tri].v0 = 1; triangles[tri].v1 = 5; triangles[tri].v2 = 4; tri++;

		// top side
		face_colors[tri] = float3(1.0f);  triangles[tri].v0 = 2; triangles[tri].v1 = 6; triangles[tri].v2 = 3; tri++;
		face_colors[tri] = float3(1.0f);  triangles[tri].v0 = 3; triangles[tri].v1 = 6; triangles[tri].v2 = 7; tri++;

		// front side
		face_colors[tri] = float3(0, 0, 1); triangles[tri].v0 = 0; triangles[tri].v1 = 4; triangles[tri].v2 = 2; tri++;
		face_colors[tri] = float3(0, 0, 1); triangles[tri].v0 = 2; triangles[tri].v1 = 4; triangles[tri].v2 = 6; tri++;

		// back side
		face_colors[tri] = float3(1, 1, 0); triangles[tri].v0 = 1; triangles[tri].v1 = 3; triangles[tri].v2 = 5; tri++;
		face_colors[tri] = float3(1, 1, 0); triangles[tri].v0 = 3; triangles[tri].v1 = 7; triangles[tri].v2 = 5; tri++;

		rtcUnmapBuffer(scene_i, mesh, RTC_INDEX_BUFFER);

		rtcSetBuffer(scene_i, mesh, RTC_USER_VERTEX_BUFFER0, vertex_colors, 0, sizeof(float3));

		return mesh;
	}

	/* adds a ground plane to the scene */
	unsigned int addGroundPlane(RTCScene scene_i)
	{
		/* create a triangulated plane with 2 triangles and 4 vertices */
		unsigned int mesh = rtcNewTriangleMesh(scene_i, RTC_GEOMETRY_STATIC, 2, 4);

		/* set vertices */
		Vertex* vertices = (Vertex*)rtcMapBuffer(scene_i, mesh, RTC_VERTEX_BUFFER);
		vertices[0].x = -10; vertices[0].y = -2; vertices[0].z = -10;
		vertices[1].x = -10; vertices[1].y = -2; vertices[1].z = +10;
		vertices[2].x = +10; vertices[2].y = -2; vertices[2].z = -10;
		vertices[3].x = +10; vertices[3].y = -2; vertices[3].z = +10;
		rtcUnmapBuffer(scene_i, mesh, RTC_VERTEX_BUFFER);

		/* set triangles */
		Triangle* triangles = (Triangle*)rtcMapBuffer(scene_i, mesh, RTC_INDEX_BUFFER);
		triangles[0].v0 = 0; triangles[0].v1 = 2; triangles[0].v2 = 1;
		triangles[1].v0 = 1; triangles[1].v1 = 2; triangles[1].v2 = 3;
		rtcUnmapBuffer(scene_i, mesh, RTC_INDEX_BUFFER);

		return mesh;
	}

	float3 castRay(float3 rayOrg, float3 rayDir)
	{
		/* initialize ray */
		RTCRay ray;
		memcpy(ray.org, begin(rayOrg), sizeof(float3));
		memcpy(ray.dir, begin(rayDir), sizeof(float3));
		ray.tnear = 0.0f;
		ray.tfar = INFINITY;
		ray.geomID = RTC_INVALID_GEOMETRY_ID;
		ray.primID = RTC_INVALID_GEOMETRY_ID;
		ray.mask = -1;
		ray.time = 0;

		/* intersect ray with scene */
		rtcIntersect(g_scene, ray);

		/* shade pixels */
		float3 color = float3(0.0f);
		if (ray.geomID != RTC_INVALID_GEOMETRY_ID)
		{
			float3 diffuse = face_colors[ray.primID];
			color = color + diffuse*0.5f;
			float3 lightDir = normalize(float3(-1, -1, -1));

			/* initialize shadow ray */
			auto shadowOrg = (rayOrg + rayDir*ray.tfar);
			auto shadowDir = (lightDir*-1.0f);
			RTCRay shadow;
			memcpy(shadow.org, begin(shadowOrg), sizeof(float3));
			memcpy(shadow.dir, begin(shadowDir), sizeof(float3));
			shadow.tnear = 0.001f;
			shadow.tfar = INFINITY;
			shadow.geomID = 1;
			shadow.primID = 0;
			shadow.mask = -1;
			shadow.time = 0;

			/* trace shadow ray */
			rtcOccluded(g_scene, shadow);

			/* add light contribution */
			if (shadow.geomID)
				color = color + diffuse*clamp(-dot(lightDir, normalize(float3(ray.Ng))), 0.0f, 1.0f);
		}
		return color;
	}

} // namespace embree

using namespace embree;
using namespace std;

void ViewImageFile(const wstring& fileName);

void main(int argc, char* argv[]) {
	// Image
	auto width = 512, height = 512;
	auto pixels = new Color[width*height];

	// Embree device and scene
	g_device = rtcNewDevice(nullptr);
	g_scene = rtcDeviceNewScene(g_device, RTC_SCENE_STATIC, RTC_INTERSECT1);
	addCube(g_scene);
	addGroundPlane(g_scene);
	rtcCommit(g_scene);

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

	// Render scene
	auto startTime = chrono::high_resolution_clock::now();
	parallel_for(size_t(0), size_t(width*height), [&](size_t i) {
		auto x = (i%width);
		auto y = (i/width);
		auto u = (2 * (x * dx) - 1);
		auto v = (2 * (y * dy) - 1);
		auto d = mul(viewMatrixInverse, float3(u * aspect, -v, z)) - from;
		auto c = castRay(from, normalize(d));
		pixels[i] = color(c);
	});
	auto endTime = chrono::high_resolution_clock::now();
	cout << "Render time = " << chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() << endl;

	// Free device and scene
	rtcDeleteScene(g_scene); g_scene = nullptr;
	rtcDeleteDevice(g_device); g_device = nullptr;
	alignedFree(face_colors); face_colors = nullptr;
	alignedFree(vertex_colors); vertex_colors = nullptr;

	// Write image file
	auto file = ofstream(L"image.ppm");
	file << "P3" << endl;
	file << width << " " << height << endl;
	file << 255 << endl;
	for (int i = 0; i < width*height; i++) {
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

void ViewImageFile(const wstring& fileName) {
	ShellExecute(0, 0, fileName.c_str(), 0, 0, SW_SHOW);
}