
**Status:** In development (already done, just getting all the files together and writing some documentation).

*Information is chronologically added to this README.md file. When enough is implemented
and fleshed out, then the library will be released.*


# LinearAlgebra

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


*[LinearAlgebra][LinearAlgebra] was first created in 2013 for use with OpenGL in Java. Later
it was rewritten for C++ and Python, where the C++ version now is the main version.*


### License

This module is shared under the MIT license, and is therefore free to use, share, distribute and modify.
See [LICENSE][LinearAlgebraLicense] for more details.


[LinearAlgebra]: https://github.com/MrVallentin/LinearAlgebra

[LinearAlgebraIssueTracker]: https://github.com/MrVallentin/LinearAlgebra/issues
[LinearAlgebraLicense]: https://github.com/MrVallentin/LinearAlgebra/blob/master/LICENSE

[LinearAlgebraWiki]: https://en.wikipedia.org/wiki/Linear_algebra
[LinearMapWiki]: https://en.wikipedia.org/wiki/Linear_map

[LinearInterpolationWiki]: https://en.wikipedia.org/wiki/Linear_interpolation

[GLSL]: https://www.opengl.org/documentation/glsl/
