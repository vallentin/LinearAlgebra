
**Status:** In development (already done, just getting all the files together and writing some documentation).

*Information is chronologically added to this README.md file. When enough is implemented
and fleshed out, then the library will be released.*


# LinearMath [![Build Status](https://drone.io/github.com/MrVallentin/LinearMath/status.png)](https://drone.io/github.com/MrVallentin/LinearMath/latest)

[LinearMath][LinearMath] is a linear algebra library, which consists of:

- Vectors
 - Vector 2D (`vec2`)
 - Vector 3D (`vec3`)
 - Vector 4D (`vec4`)

- Matrices
 - Matrix 3D (`mat3`, `mat3x3`)
 - Matrix 4D (`mat4`, `mat4x4`)

- Quaternion (`quat`)
- Dual Quaternion (`dualquat`)

Since [LinearMath][LinearMath] was built for use with computer graphics, it contains
a few helper classes for especially that:

- Transform
- Project
- MatrixStack

*[LinearMath][LinearMath] was first created in 2013 for use with OpenGL in Java. Later
it was remade for C++ and Python, where the C++ version now is the main version.*


### License

This module is shared under the MIT license, and is therefore free to use, share, distribute and modify.
See [LICENSE][LinearMathLicense] for more details.


[LinearMath]: https://github.com/MrVallentin/LinearMath

[LinearMathIssueTracker]: https://github.com/MrVallentin/LinearMath/issues
[LinearMathLicense]: https://github.com/MrVallentin/LinearMath/blob/master/LICENSE

[LinearAlgebra]: https://en.wikipedia.org/wiki/Linear_algebra
[LinearMap]: https://en.wikipedia.org/wiki/Linear_map

[LinearInterpolation]: https://en.wikipedia.org/wiki/Linear_interpolation

[GLSL]: https://www.opengl.org/documentation/glsl/
