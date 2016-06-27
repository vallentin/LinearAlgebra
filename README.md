
# LinearAlgebra [![Build Status][LinearAlgebraBuildStatus]][LinearAlgebraCI] ![Release][LinearAlgebraVersionBadge] ![License][LinearAlgebraLicenseBadge]

[LinearAlgebra][LinearAlgebra] is a linear algebra library, which consists of:

- Vectors
 - Vector 2D (`vec2`)
 - Vector 3D (`vec3`)
 - Vector 4D (`vec4`)

- Matrices
 - Matrix 2D (`mat2`, `mat2x2`)
 - Matrix 3D (`mat3`, `mat3x3`)
 - Matrix 4D (`mat4`, `mat4x4`)

- Quaternion (`quat`)
- Dual Quaternion (`dualquat`)


Since [LinearAlgebra][LinearAlgebra] was built for use with computer graphics, it contains
a few helper classes for especially that:

- Transform
- Project
- MatrixStack


Various methods are implemented such that they can be called by your liking:

```cpp
vec3 a(1.0f, 2.0f, 3.0f);
vec3 b(4.0f, 5.0f, 6.0f);

float dotProduct = a.dot(b);
// vs
float dotProduct = dot(a, b);
```

```cpp
mat4 model = mat4::identity;
// ...

mat3 normalMatrix = model.inverse().transpose();
// vs
mat3 normalMatrix = transpose(inverse(model));
```


If GCC is used then compiling [LinearAlgebra][LinearAlgebra] must include `-Wno-unknown-pragmas`,
as various `#pragma region`'s are scattered in [LinearAlgebra][LinearAlgebra]'s implementation.


*[LinearAlgebra][LinearAlgebra] was first created in 2013 for use with OpenGL in Java. Later
it was rewritten for C++ and Python, where the C++ version now is the main version.*


### License

This module is shared under the MIT license, and is therefore free to use, share, distribute and modify.
See [LICENSE][LinearAlgebraLicense] for more details.


[LinearAlgebra]: https://github.com/MrVallentin/LinearAlgebra
[LinearAlgebraLicense]: https://github.com/MrVallentin/LinearAlgebra/blob/master/LICENSE

[LinearAlgebraBuildStatus]: https://drone.io/github.com/MrVallentin/LinearAlgebra/status.png
[LinearAlgebraCI]: https://drone.io/github.com/MrVallentin/LinearAlgebra/latest

[LinearAlgebraVersionBadge]: https://img.shields.io/badge/release-v0.1.0-blue.svg
[LinearAlgebraLicenseBadge]: https://img.shields.io/badge/license-%20free%20to%20use%2C%20share%2C%20modify%20and%20redistribute-blue.svg

[LinearAlgebraIssueTracker]: https://github.com/MrVallentin/LinearAlgebra/issues


[LinearAlgebraWiki]: https://en.wikipedia.org/wiki/Linear_algebra
[LinearMapWiki]: https://en.wikipedia.org/wiki/Linear_map

[LinearInterpolationWiki]: https://en.wikipedia.org/wiki/Linear_interpolation

[GLSL]: https://www.opengl.org/documentation/glsl/
