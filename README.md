# EmbreeTest
This is a simple standalone version of the [Intel Embree](https://github.com/embree/embree) [triangle_geometry tutorial](https://github.com/embree/embree/blob/master/tutorials/triangle_geometry/triangle_geometry_device.cpp).

![](Result.png?raw=true)

Install the [Intel Embree v2.15.0 x64 Windows binary distribution](https://embree.github.io/downloads.html).

Build with Visual Studio 2017.

If you install another version of Embree, or put it somewhere other than your C:\ drive, you'll need to edit **Embree.props** and change the **EMBREEx64** path to be the new location. 

The generated image is in PPM format, inspired by [Kevin Beason's smallpt 99 line renderer](http://www.kevinbeason.com/smallpt/). The viewer I use for PPM files is [DJV](http://djv.sourceforge.net/index.html).

I created this sample as I found it painful to compile the whole Embree source tree just to try out the tutorials.

Embree uses its own internal vector and matrix library, which aren't included when you install the binary distribution, so I used the single header file template library [linalg.h](https://github.com/sgorsten/linalg) instead. 

This code is released under the license is Apache 2.0, the same as the original Intel license.
