
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
- ~~Dual Quaternion (`dualquat`)~~


Since [LinearAlgebra][LinearAlgebra] was built for use with computer graphics, it contains
a few helper classes for especially that:

- ~~Transform~~
- ~~Project~~ (Functionality stored in `mat4`)
- ~~MatrixStack~~


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


### Mouse to Ray

```cpp
vec2 mousePos = ...;

ivec4 viewport = ivec4(0, 0, viewportWidth, viewportHeight);

mat4 view = ...;
mat4 projection = ...;

vec3 rayStart = mat4::unproject(vec3(mousePos, 0.0f), view, projection, viewport);
vec3 rayEnd   = mat4::unproject(vec3(mousePos, 1.0f), view, projection, viewport);
vec3 rayDirection = normalize(rayEnd - rayStart);
```


### std::cout & std::cin

To add friend functions for stream purposes,
simply make sure to have the `iostream` header
included prior to including `linalg`.

```cpp
#include <iostream>
#include "linalg.hpp"
```

This now allows printing and reading `vec`, `mat` and `quat`.

```cpp
std::cout << vec4(1.0f, 2.0f, 3.0f, 4.0f) << std::endl;

// Prints:
// vec4(1.0, 2.0, 3.0, 4.0)

std::cout << mat3(
	1.0f, 2.0f, 3.0f,
	4.0f, 5.0f, 6.0f,
	7.0f, 8.0f, 9.0f) << std::endl;

// Prints:
// mat3(1.0, 2.0, 3.0
//      4.0, 5.0, 6.0,
//      7.0, 8.0, 9.0)
```


## License

This module is shared under the MIT license, and is therefore free to use, share, distribute and modify.
See [LICENSE][LinearAlgebraLicense] for more details.


[LinearAlgebra]: https://github.com/MrVallentin/LinearAlgebra
[LinearAlgebraLicense]: https://github.com/MrVallentin/LinearAlgebra/blob/master/LICENSE

[LinearAlgebraBuildStatus]: https://drone.io/github.com/MrVallentin/LinearAlgebra/status.png
[LinearAlgebraCI]: https://drone.io/github.com/MrVallentin/LinearAlgebra/latest

[LinearAlgebraVersionBadge]: https://img.shields.io/badge/release-v1.1.17-blue.svg
[LinearAlgebraLicenseBadge]: https://img.shields.io/badge/license-%20free%20to%20use%2C%20share%2C%20modify%20and%20redistribute-blue.svg

[LinearAlgebraIssueTracker]: https://github.com/MrVallentin/LinearAlgebra/issues


[LinearAlgebraWiki]: https://en.wikipedia.org/wiki/Linear_algebra
[LinearMapWiki]: https://en.wikipedia.org/wiki/Linear_map

[LinearInterpolationWiki]: https://en.wikipedia.org/wiki/Linear_interpolation

[GLSL]: https://www.opengl.org/documentation/glsl/
