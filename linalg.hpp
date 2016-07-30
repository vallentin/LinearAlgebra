
// Author: Christian Vallentin <mail@vallentinsource.com>
// Website: http://vallentinsource.com
// Repository: https://github.com/MrVallentin/LinearAlgebra
//
// Date Created: October 01, 2013
// Last Modified: July 16, 2016

// Copyright (c) 2013-2016 Christian Vallentin <mail@vallentinsource.com>
//
// This software is provided 'as-is', without any express or implied
// warranty. In no event will the authors be held liable for any damages
// arising from the use of this software.
//
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
//
// 1. The origin of this software must not be misrepresented; you must not
//    claim that you wrote the original software. If you use this software
//    in a product, an acknowledgment in the product documentation would be
//    appreciated but is not required.
//
// 2. Altered source versions must be plainly marked as such, and must not
//    be misrepresented as being the original software.
//
// 3. This notice may not be removed or altered from any source
//    distribution.

// Refrain from using any exposed macros, functions
// or structs prefixed with an underscore. As these
// are only intended for internal purposes. Which
// additionally means they can be removed, renamed
// or changed between minor updates without notice.

#ifndef LINEAR_ALGEBRA_HPP
#define LINEAR_ALGEBRA_HPP


#define _LINALG_STRINGIFY(str) #str
#define _LINALG_STRINGIFY_TOKEN(str) _LINALG_STRINGIFY(str)

#define LINALG_STRINGIFY_VERSION(major, minor, patch) _LINALG_STRINGIFY(major) "." _LINALG_STRINGIFY(minor) "." _LINALG_STRINGIFY(patch)


#define LINALG_NAME "LinearAlgebra"

#define LINALG_VERSION_MAJOR 1
#define LINALG_VERSION_MINOR 1
#define LINALG_VERSION_PATCH 17

#define LINALG_VERSION LINALG_STRINGIFY_VERSION(LINALG_VERSION_MAJOR, LINALG_VERSION_MINOR, LINALG_VERSION_PATCH)

#define LINALG_NAME_VERSION LINALG_NAME " " LINALG_VERSION


#include <math.h>


#ifdef _IOSTREAM_

#define _LINALG_IN_FIND_BEGIN() \
	do { \
		stream >> c; \
		if (stream.eof()) \
			return stream; \
	} while (c != '(');

#define _LINALG_IN_FIND_NEXT() \
	do { \
		stream >> c; \
		if (stream.eof()) \
			return stream; \
		else if (c == ')') \
			return stream; \
	} while (c != ',');

#define _LINALG_IN_FIND_END() \
	do { \
		stream >> c; \
		if (stream.eof()) \
			return stream; \
	} while (c != ')');

#endif


#ifndef LINALG_DEFAULT_SCALAR
#	define LINALG_DEFAULT_SCALAR float
#endif


// Disable structure padding
#pragma pack(push, 1)


template<typename T> class vec2_t;
template<typename T> class vec3_t;
template<typename T> class vec4_t;

template<typename T> class mat2_t;
template<typename T> class mat3_t;
template<typename T> class mat4_t;

template<typename T> class quat_t;


typedef vec2_t<LINALG_DEFAULT_SCALAR> vec2;

typedef vec2_t<float> fvec2;
typedef vec2_t<double> dvec2;

typedef vec2_t<signed int> ivec2;
typedef vec2_t<unsigned int> uvec2;

typedef vec2_t<bool> bvec2;

typedef vec2_t<signed long> lvec2;
typedef vec2_t<unsigned long> ulvec2;

typedef vec2_t<signed long long> llvec2;
typedef vec2_t<unsigned long long> ullvec2;


typedef vec3_t<LINALG_DEFAULT_SCALAR> vec3;

typedef vec3_t<float> fvec3;
typedef vec3_t<double> dvec3;

typedef vec3_t<signed int> ivec3;
typedef vec3_t<unsigned int> uvec3;

typedef vec3_t<bool> bvec3;

typedef vec3_t<signed long> lvec3;
typedef vec3_t<unsigned long> ulvec3;

typedef vec3_t<signed long long> llvec3;
typedef vec3_t<unsigned long long> ullvec3;


typedef vec4_t<LINALG_DEFAULT_SCALAR> vec4;

typedef vec4_t<float> fvec4;
typedef vec4_t<double> dvec4;

typedef vec4_t<signed int> ivec4;
typedef vec4_t<unsigned int> uvec4;

typedef vec4_t<bool> bvec4;

typedef vec4_t<signed long> lvec4;
typedef vec4_t<unsigned long> ulvec4;

typedef vec4_t<signed long long> llvec4;
typedef vec4_t<unsigned long long> ullvec4;


typedef mat2_t<LINALG_DEFAULT_SCALAR> mat2;

typedef mat2_t<float> fmat2;
typedef mat2_t<double> dmat2;

template<typename T> using mat2x2_t = mat2_t<T>;

typedef mat2x2_t<LINALG_DEFAULT_SCALAR> mat2x2;

typedef mat2x2_t<float> fmat2x2;
typedef mat2x2_t<double> dmat2x2;


typedef mat3_t<LINALG_DEFAULT_SCALAR> mat3;

typedef mat3_t<float> fmat3;
typedef mat3_t<double> dmat3;

template<typename T> using mat3x3_t = mat3_t<T>;

typedef mat3x3_t<LINALG_DEFAULT_SCALAR> mat3x3;

typedef mat3x3_t<float> fmat3x3;
typedef mat3x3_t<double> dmat3x3;


typedef mat4_t<LINALG_DEFAULT_SCALAR> mat4;

typedef mat4_t<float> fmat4;
typedef mat4_t<double> dmat4;

template<typename T> using mat4x4_t = mat4_t<T>;

typedef mat4x4_t<LINALG_DEFAULT_SCALAR> mat4x4;

typedef mat4x4_t<float> fmat4x4;
typedef mat4x4_t<double> dmat4x4;


typedef quat_t<LINALG_DEFAULT_SCALAR> quat;

typedef quat_t<float> fquat;
typedef quat_t<double> dquat;


#if defined(_DEBUG) && !defined(DEBUG)
#	define DEBUG 1
#endif


// This was changed from 1E-6 to 1E-4, as asserting rotate(90deg) didn't match.
#define LINALG_EPSILON 1E-4f

#define LINALG_FEQUAL(x, y) ((((y) - LINALG_EPSILON) < (x)) && ((x) < ((y) + LINALG_EPSILON)))
#define LINALG_DEQUAL(x, y) ((((y) - LINALG_EPSILON) < (x)) && ((x) < ((y) + LINALG_EPSILON)))

#define LINALG_PI 3.1415926535897932

#define LINALG_DEG2RAD (LINALG_PI / 180.0)
#define LINALG_RAD2DEG (180.0 / LINALG_PI)


// These annoyingly named Windows macros are interfering with the related vector methods!

#ifdef min
#	undef min
#endif

#ifdef max
#	undef max
#endif


template<typename T>
class vec2_t
{
private:

	typedef vec2_t<T> vec2;
	typedef vec3_t<T> vec3;
	typedef vec4_t<T> vec4;


public:

	static const vec2_t<T> zero;
	static const vec2_t<T> one;

	static const vec2_t<T> up, down;
	static const vec2_t<T> left, right;


public:

	// Performs Gram-Schmidt Orthogonalization on 2 basis vectors to turn them into orthonormal basis vectors
	static vec2 orthogonalize(const vec2 &a, vec2 &b)
	{
		b = b - b.project(a);
		b = b.normalize();
	}


public:

	T x, y;


public:

	vec2_t() : x(T(0)), y(T(0)) {}

	vec2_t(const vec2_t<T> &v) : x(T(v.x)), y(T(v.y)) {}
	template<typename T2> vec2_t(const vec2_t<T2> &v) : x(T(v.x)), y(T(v.y)) {}

	template<typename T2> vec2_t(const T2 &xy) : x(T(xy)), y(T(xy)) {}
	template<typename T2> vec2_t(const T2 &x, const T2 &y) : x(T(x)), y(T(y)) {}

	template<typename T2> vec2_t(const T2 *xy) : x(T(xy[0])), y(T(xy[1])) {}

	template<typename T2> vec2_t(const vec3_t<T2> &v);
	template<typename T2> vec2_t(const vec4_t<T2> &v);

	~vec2_t() {}


#pragma region Operator Overloading

#pragma region Member Access Operators

	inline T& operator[](const int index) { return (reinterpret_cast<T*>(this))[index]; }
	inline T operator[](const int index) const { return ((T*) this)[index]; }

#pragma endregion

#pragma region Arithmetic Operators

	vec2 operator+() const { return vec2(+this->x, +this->y); }
	vec2 operator-() const { return vec2(-this->x, -this->y); }

	friend vec2 operator+(const vec2 &lhs, const vec2 &rhs) { return vec2((lhs.x + rhs.x), (lhs.y + rhs.y)); }
	friend vec2 operator-(const vec2 &lhs, const vec2 &rhs) { return vec2((lhs.x - rhs.x), (lhs.y - rhs.y)); }
	friend vec2 operator*(const vec2 &lhs, const vec2 &rhs) { return vec2((lhs.x * rhs.x), (lhs.y * rhs.y)); }
	friend vec2 operator/(const vec2 &lhs, const vec2 &rhs) { return vec2((lhs.x / rhs.x), (lhs.y / rhs.y)); }
	friend vec2 operator%(const vec2 &lhs, const vec2 &rhs) { return vec2((lhs.x % rhs.x), (lhs.y % rhs.y)); }

	friend inline vec2 operator+(const vec2 &lhs, const T &rhs) { return (lhs + vec2(rhs)); }
	friend inline vec2 operator+(const T &lhs, const vec2 &rhs) { return (vec2(lhs) + rhs); }

	friend inline vec2 operator-(const vec2 &lhs, const T &rhs) { return (lhs - vec2(rhs)); }
	friend inline vec2 operator-(const T &lhs, const vec2 &rhs) { return (vec2(lhs) - rhs); }

	friend inline vec2 operator*(const vec2 &lhs, const T &rhs) { return (lhs * vec2(rhs)); }
	friend inline vec2 operator*(const T &lhs, const vec2 &rhs) { return (vec2(lhs) * rhs); }

	friend inline vec2 operator/(const vec2 &lhs, const T &rhs) { return (lhs / vec2(rhs)); }
	friend inline vec2 operator/(const T &lhs, const vec2 &rhs) { return (vec2(lhs) / rhs); }

	friend inline vec2 operator%(const vec2 &lhs, const T &rhs) { return (lhs % vec2(rhs)); }
	friend inline vec2 operator%(const T &lhs, const vec2 &rhs) { return (vec2(lhs) % rhs); }

#pragma endregion
#pragma region Increment & Decrement Operators

	vec2& operator++() // Prefix
	{
		++this->x;
		++this->y;

		return (*this);
	}
	inline vec2 operator++(int) { return ++(*this); } // Postfix

	vec2& operator--() // Prefix
	{
		--this->x;
		--this->y;

		return (*this);
	}
	inline vec2 operator--(int) { return --(*this); } // Postfix

#pragma endregion
#pragma region Assignment Operators

	inline vec2& operator+=(const vec2 &rhs) { return ((*this) = ((*this) + rhs)); }
	inline vec2& operator+=(const T &rhs) { return ((*this) = ((*this) + rhs)); }

	inline vec2& operator-=(const vec2 &rhs) { return ((*this) = ((*this) - rhs)); }
	inline vec2& operator-=(const T &rhs) { return ((*this) = ((*this) - rhs)); }

	inline vec2& operator*=(const vec2 &rhs) { return ((*this) = ((*this) * rhs)); }
	inline vec2& operator*=(const T &rhs) { return ((*this) = ((*this) * rhs)); }

	inline vec2& operator/=(const vec2 &rhs) { return ((*this) = ((*this) / rhs)); }
	inline vec2& operator/=(const T &rhs) { return ((*this) = ((*this) / rhs)); }

	inline vec2& operator%=(const vec2 &rhs) { return ((*this) = ((*this) % rhs)); }
	inline vec2& operator%=(const T &rhs) { return ((*this) = ((*this) % rhs)); }

	vec2& operator=(const T &rhs)
	{
		(*this) = vec2(rhs);

		return (*this);
	}

	template<typename T2>
	vec2& operator=(const T2 &rhs)
	{
		(*this) = vec2(rhs);

		return (*this);
	}

	vec2& operator=(const vec2 &rhs)
	{
		this->x = rhs.x;
		this->y = rhs.y;

		return (*this);
	}

	template<typename T2>
	vec2& operator=(const vec2_t<T2> &rhs)
	{
		this->x = T(rhs.x);
		this->y = T(rhs.y);

		return (*this);
	}

#pragma endregion

#pragma region Logical Operators
#pragma endregion
#pragma region Comparison Operators

	bool operator==(const vec2 &rhs) const;

	friend inline bool operator==(const vec2 &lhs, const T &rhs) { return (lhs == vec2(rhs)); }
	friend inline bool operator==(const T &lhs, const vec2 &rhs) { return (vec2(lhs) == rhs); }

	friend inline bool operator!=(const vec2 &lhs, const vec2 &rhs) { return !(lhs == rhs); }
	friend inline bool operator!=(const vec2 &lhs, const T &rhs) { return (lhs != vec2(rhs)); }
	friend inline bool operator!=(const T &lhs, const vec2 &rhs) { return (vec2(lhs) != rhs); }

	friend inline bool operator>(const vec2 &lhs, const vec2 &rhs) { return ((lhs.x > rhs.x) && (lhs.y > rhs.y)); }
	friend inline bool operator>=(const vec2 &lhs, const vec2 &rhs) { return ((lhs.x >= rhs.x) && (lhs.y >= rhs.y)); }
	friend inline bool operator<(const vec2 &lhs, const vec2 &rhs) { return ((lhs.x < rhs.x) && (lhs.y < rhs.y)); }
	friend inline bool operator<=(const vec2 &lhs, const vec2 &rhs) { return ((lhs.x <= rhs.x) && (lhs.y <= rhs.y)); }

	friend inline bool operator>(const vec2 &lhs, const T &rhs) { return (lhs > vec2(rhs)); }
	friend inline bool operator>(const T &lhs, const vec2 &rhs) { return (vec2(lhs) > rhs); }

	friend inline bool operator>=(const vec2 &lhs, const T &rhs) { return (lhs >= vec2(rhs)); }
	friend inline bool operator>=(const T &lhs, const vec2 &rhs) { return (vec2(lhs) >= rhs); }

	friend inline bool operator<(const vec2 &lhs, const T &rhs) { return (lhs < vec2(rhs)); }
	friend inline bool operator<(const T &lhs, const vec2 &rhs) { return (vec2(lhs) < rhs); }

	friend inline bool operator<=(const vec2 &lhs, const T &rhs) { return (lhs <= vec2(rhs)); }
	friend inline bool operator<=(const T &lhs, const vec2 &rhs) { return (vec2(lhs) <= rhs); }

#pragma endregion

#pragma region Cast Operators

	inline explicit operator T*() const
	{
		return reinterpret_cast<T*>(this);
	}

	template<typename T2>
	inline operator vec2_t<T2>() const
	{
		return vec2_t<T2>(
			static_cast<T2>(this->x),
			static_cast<T2>(this->y)
		);
	}

#pragma endregion

#pragma region Stream Operators

#ifdef _IOSTREAM_

	friend inline std::ostream& operator<<(std::ostream &stream, const vec2 &rhs)
	{
		return (stream << "vec2(" << rhs.x << ", " << rhs.y << ")");
	}

	friend inline std::wostream& operator<<(std::wostream &stream, const vec2 &rhs)
	{
		return (stream << L"vec2(" << rhs.x << L", " << rhs.y << L")");
	}

	friend inline std::istream& operator>>(std::istream &stream, vec2 &rhs)
	{
		rhs = vec2::zero;

		char c;

		_LINALG_IN_FIND_BEGIN();
		stream >> rhs.x;

		_LINALG_IN_FIND_NEXT();
		stream >> rhs.y;

		_LINALG_IN_FIND_END();

		return stream;
	}

	friend inline std::wistream& operator>>(std::wistream &stream, vec2 &rhs)
	{
		rhs = vec2::zero;

		wchar_t c;

		_LINALG_IN_FIND_BEGIN();
		stream >> rhs.x;

		_LINALG_IN_FIND_NEXT();
		stream >> rhs.y;

		_LINALG_IN_FIND_END();

		return stream;
	}

#endif

#pragma endregion

#pragma endregion

#pragma region

	inline T dot(const vec2 &rhs) const
	{
		return (this->x * rhs.x + this->y * rhs.y);
	}
	friend inline T dot(const vec2 &lhs, const vec2 &rhs) { return lhs.dot(rhs); }


	vec2 cross(const vec2 &rhs) const
	{
		return vec2(
			((this->y * rhs.z) - (this->z * rhs.y)),
			((this->z * rhs.x) - (this->x * rhs.z))
		);
	}
	friend inline vec2 cross(const vec2 &lhs, const vec2 &rhs) { return lhs.cross(rhs); }


	inline T lengthSquared() const { return (this->x * this->x + this->y * this->y); }
	friend inline T lengthSquared(const vec2 &lhs) { return lhs.lengthSquared(); }

	inline T length2() const { return (this->x * this->x + this->y * this->y); }
	friend inline T length2(const vec2 &lhs) { return lhs.lengthSquared(); }

	inline T length() const { return sqrt(lengthSquared()); }
	friend inline T length(const vec2 &lhs) { return lhs.length(); }


	inline T distanceSquared(const vec2 &rhs) { return ((*this) - rhs).lengthSquared(); }
	friend inline T distanceSquared(const vec2 &lhs, const vec2 &rhs) { return lhs.distanceSquared(rhs); }

	inline T distance(const vec2 &rhs) { return ((*this) - rhs).length(); }
	friend inline T distance(const vec2 &lhs, const vec2 &rhs) { return lhs.distance(rhs); }


	vec2 normalize(const T &to = T(1.0)) const;
	friend inline vec2 normalize(const vec2 &lhs, const T &to = T(1.0)) { return lhs.normalize(to); }


	T angle(const vec2 &rhs) const
	{
		return (dot(rhs) / (length() * rhs.length()));
	}
	friend inline T angle(const vec2 &lhs, const vec2 &rhs) { return lhs.angle(rhs); }


	// Calculates the projection of a onto b
	// 
	// Reference: http://en.wikipedia.org/wiki/Vector_projection#Vector_projection_2
	inline vec2 project(const vec2 &b) const
	{
		const float length = b.length();

		return ((dot(b) / (length * length)) * b);
	}
	friend inline vec2 project(const vec2 &a, const vec2 &b) { return a.project(b); }


	// Calculates the components of a perpendicular to b
	inline vec2 perpendicular(const vec2 &b) const
	{
		const float length = b.length();

		return ((*this) - ((dot(b) / (length * length)) * b));
	}
	friend inline vec2 perpendicular(const vec2 &a, const vec2 &b) { return a.perpendicular(b); }


	// Calculates the reflection vector from entering ray direction a and surface normal b
	inline vec2 reflect(const vec2 &b) const
	{
		return (T(2) * project(b) - (*this));
	}
	friend inline vec2 reflect(const vec2 &a, const vec2 &b) { return a.reflect(b); }


	inline T cosine(const vec2 &b) const
	{
		return normalize().dot(b.normalize());
	}
	friend inline T cosine(const vec2 &a, const vec2 &b) { return a.cosine(b); }


	inline vec2 rotate(const T theta) const
	{
		const T c = cos(theta), s = sin(theta);

		return vec2(
			(c * this->x - s * this->y),
			(s * this->x + c * this->y)
		);
	}
	friend inline vec2 rotate(const vec2 &a, const T theta) { return a.rotate(theta); }

#pragma endregion

#pragma region

	inline bool isNullVector() const;
	friend inline bool isNullVector(const vec2 &v) { return v.isNullVector(); }

	inline bool isUnitVector() const;
	friend inline bool isUnitVector(const vec2 &v) { return v.isUnitVector(); }

	inline bool isNormalized(const T &to = T(1.0)) const;
	friend inline bool isNormalized(const vec2 &v) { return v.isNormalized(); }


	inline bool isOrthogonalTo(const vec2 &rhs) const;
	friend inline bool isOrthogonalTo(const vec2 &lhs, const vec2 &rhs) { return (lhs.isOrthogonalTo(rhs)); }

	inline bool isPerpendicularTo(const vec2 &rhs) const;
	friend inline bool isPerpendicularTo(const vec2 &lhs, const vec2 &rhs) { return (lhs.isPerpendicularTo(rhs)); }

	inline bool isParallelTo(const vec2 &rhs) const;
	friend inline bool isParallelTo(const vec2 &lhs, const vec2 &rhs) { return (lhs.isParallelTo(rhs)); }

#pragma endregion

#pragma region

	inline vec2 abs() const { return vec2(abs(this->x), abs(this->y), abs(this->z)); }
	friend inline vec2 abs(const vec2 &v) { return v.abs(); }


	inline vec2 max(const vec2 &rhs) const
	{
		return vec2(
			((this->x > rhs.x) ? this->x : rhs.x),
			((this->y > rhs.y) ? this->y : rhs.y)
		);
	}
	friend inline vec2 max(const vec2 &lhs, const vec2 &rhs) { return lhs.max(rhs); }

	inline vec2 min(const vec2 &rhs) const
	{
		return vec2(
			((this->x < rhs.x) ? this->x : rhs.x),
			((this->y < rhs.y) ? this->y : rhs.y)
		);
	}
	friend inline vec2 min(const vec2 &lhs, const vec2 &rhs) { return lhs.min(rhs); }

	inline vec2 clamp(const vec2 &min, const vec2 &max) const
	{
		return min(max).max(min);
	}
	friend inline vec2 clamp(const vec2 &v, const vec2 &min, const vec2 &max) { return v.clamp(min, max); }


	inline vec2 lerp(const vec2 &to, const T &t) const { return ((*this) + t * (to - (*this))); }
	friend inline vec2 lerp(const vec2 &from, const vec2 &to, const T &t) { return from.lerp(to, t); }

	inline vec2 lerp(const vec2 &to, const vec2 &t) const { return ((*this) + t * (to - (*this))); }
	friend inline vec2 lerp(const vec2 &from, const vec2 &to, const vec2 &t) { return from.lerp(to, t); }


	// Reference: https://en.wikipedia.org/wiki/Slerp
	// Referemce: https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
	inline vec2 slerp(const vec2 &to, const T &t) const
	{
		const T dot = normalize().dot(to.normalize());
		const T theta = acos(dot);
		const T s = sin(theta);

		return (sin((T(1) - t) * theta) / s * (*this) + sin(t * theta) / s * to);
	}
	friend inline vec2 slerp(const vec2 &from, const vec2 &to, const T &t) { return from.slerp(to, t); }

	// Reference: https://en.wikipedia.org/wiki/Slerp
	// Referemce: https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
	inline vec2 slerp(const vec2 &to, const vec2 &t) const
	{
		const T dot = normalize().dot(to.normalize());
		const T theta = acos(dot);
		const T s = sin(theta);

		return (sin((T(1) - t) * theta) / s * (*this) + sin(t * theta) / s * to);
	}
	friend inline vec2 slerp(const vec2 &from, const vec2 &to, const vec2 &t) { return from.slerp(to, t); }


	// Reference: https://en.wikipedia.org/wiki/Sign_function
	inline const vec2 signum() const
	{
		return vec2(
			((this->x < T(0)) ? T(-1) : ((this->x > T(0)) ? T(1) : T(0))),
			((this->y < T(0)) ? T(-1) : ((this->y > T(0)) ? T(1) : T(0)))
		);
	}
	friend inline vec2 signum(const vec2 &v) { return v.signum(); }

#pragma endregion

#pragma region Swizzling

#define LINALG_SWIZZLE_INDEX(index, c) \
		if ((c == 'X') || (c == 'X') || (c == 'r') || (c == 'R') || (c == 's') || (c == 'S')) index = 0; \
		else if ((c == 'y') || (c == 'Y') || (c == 'g') || (c == 'G') || (c == 't') || (c == 'T')) index = 1; \
		else index = 0;

	vec2 swizzle(const char x, const char y) const
	{
		int x_index = 0, y_index = 1;

		LINALG_SWIZZLE_INDEX(x_index, x);
		LINALG_SWIZZLE_INDEX(y_index, y);

		return vec2((*this)[x_index], (*this)[y_index]);
	}

#undef LINALG_SWIZZLE_INDEX

#pragma endregion


	inline void swap(vec2 &other)
	{
		const vec2 tmp(*this);
		(*this) = other;
		other = tmp;
	}
	friend inline void swap(vec2 &a, vec2 &b) { a.swap(b); }
};


template<typename T>
class vec3_t
{
private:

	typedef vec2_t<T> vec2;
	typedef vec3_t<T> vec3;
	typedef vec4_t<T> vec4;


public:

	static const vec3_t<T> zero;
	static const vec3_t<T> one;

	static const vec3_t<T> up, down;
	static const vec3_t<T> left, right;
	static const vec3_t<T> forward, backward;


public:

	// Performs Gram-Schmidt Orthogonalization on 2 basis vectors to turn them into orthonormal basis vectors
	static inline vec3 orthogonalize(const vec3 &a, vec3 &b)
	{
		b = b - b.project(a);
		b = b.normalize();
	}

	// Performs Gram-Schmidt Orthogonalization on 3 basis vectors to turn them into orthonormal basis vectors
	static inline void orthogonalize(const vec3 &a, vec3 &b, vec3 &c)
	{
		b = b - b.project(a);
		b = b.normalize();

		c = c - c.project(a) - c.project(b);
		c = c.normalize();
	}


	// The order of the vertices given will affect the direction of the resulting normal.
	// The front face of the triangle is considered to be the counter-clockwise order of
	// the three vertices given.
	static inline vec3 surfaceNormal(const vec3 &point1, const vec3 &point2, const vec3 &point3)
	{
		const vec3 u = point3 - point1;
		const vec3 v = point2 - point1;

		const vec3 normal = u.cross(v);

		return normal.normalize();
	}


public:

	T x, y, z;


public:

	vec3_t() : x(T(0)), y(T(0)), z(T(0)) {}

	vec3_t(const vec3_t<T> &v) : x(v.x), y(v.y), z(v.z) {}
	template<typename T2> vec3_t(const vec3_t<T2> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)) {}

	template<typename T2> vec3_t(const T2 &xyz) : x(T(xyz)), y(T(xyz)), z(T(xyz)) {}
	template<typename T2> vec3_t(const T2 &x, const T2 &y, const T2 &z = T(0)) : x(T(x)), y(T(y)), z(T(z)) {}

	template<typename T2> vec3_t(const T2 *xyz) : x(T(xyz[0])), y(T(xyz[1])), z(T(xyz[2])) {}

	template<typename T2, typename T3> vec3_t(const vec2_t<T2> &xy, const T3 &z = T3(0)) : x(T(xy.x)), y(T(xy.y)), z(T(z)) {}
	template<typename T2, typename T3> vec3_t(const T2 &x, const vec2_t<T3> &yz) : x(T(x)), y(T(yz.x)), z(T(yz.y)) {}

	template<typename T2> vec3_t(const vec4_t<T2> &v);

	~vec3_t() {}


#pragma region Operator Overloading

#pragma region Member Access Operators

	inline T& operator[](const int index) { return (reinterpret_cast<T*>(this))[index]; }
	inline T operator[](const int index) const { return ((T*) this)[index]; }

#pragma endregion

#pragma region Arithmetic Operators

	vec3 operator+() const { return vec3(+this->x, +this->y, +this->z); }
	vec3 operator-() const { return vec3(-this->x, -this->y, -this->z); }

	friend vec3 operator+(const vec3 &lhs, const vec3 &rhs) { return vec3((lhs.x + rhs.x), (lhs.y + rhs.y), (lhs.z + rhs.z)); }
	friend vec3 operator-(const vec3 &lhs, const vec3 &rhs) { return vec3((lhs.x - rhs.x), (lhs.y - rhs.y), (lhs.z - rhs.z)); }
	friend vec3 operator*(const vec3 &lhs, const vec3 &rhs) { return vec3((lhs.x * rhs.x), (lhs.y * rhs.y), (lhs.z * rhs.z)); }
	friend vec3 operator/(const vec3 &lhs, const vec3 &rhs) { return vec3((lhs.x / rhs.x), (lhs.y / rhs.y), (lhs.z / rhs.z)); }
	friend vec3 operator%(const vec3 &lhs, const vec3 &rhs) { return vec3((lhs.x % rhs.x), (lhs.y % rhs.y), (lhs.z % rhs.z)); }

	friend inline vec3 operator+(const vec3 &lhs, const T &rhs) { return (lhs + vec3(rhs)); }
	friend inline vec3 operator+(const T &lhs, const vec3 &rhs) { return (vec3(lhs) + rhs); }

	friend inline vec3 operator-(const vec3 &lhs, const T &rhs) { return (lhs - vec3(rhs)); }
	friend inline vec3 operator-(const T &lhs, const vec3 &rhs) { return (vec3(lhs) - rhs); }

	friend inline vec3 operator*(const vec3 &lhs, const T &rhs) { return (lhs * vec3(rhs)); }
	friend inline vec3 operator*(const T &lhs, const vec3 &rhs) { return (vec3(lhs) * rhs); }

	friend inline vec3 operator/(const vec3 &lhs, const T &rhs) { return (lhs / vec3(rhs)); }
	friend inline vec3 operator/(const T &lhs, const vec3 &rhs) { return (vec3(lhs) / rhs); }

	friend inline vec3 operator%(const vec3 &lhs, const T &rhs) { return (lhs % vec3(rhs)); }
	friend inline vec3 operator%(const T &lhs, const vec3 &rhs) { return (vec3(lhs) % rhs); }

#pragma endregion
#pragma region Increment & Decrement Operators

	vec3& operator++() // Prefix
	{
		++this->x;
		++this->y;
		++this->z;

		return (*this);
	}
	inline vec3 operator++(int) { return ++(*this); } // Postfix

	vec3& operator--() // Prefix
	{
		--this->x;
		--this->y;
		--this->z;

		return (*this);
	}
	inline vec3 operator--(int) { return --(*this); } // Postfix

#pragma endregion
#pragma region Assignment Operators

	inline vec3& operator+=(const vec3 &rhs) { return ((*this) = ((*this) + rhs)); }
	inline vec3& operator+=(const T &rhs) { return ((*this) = ((*this) + rhs)); }

	inline vec3& operator-=(const vec3 &rhs) { return ((*this) = ((*this) - rhs)); }
	inline vec3& operator-=(const T &rhs) { return ((*this) = ((*this) - rhs)); }

	inline vec3& operator*=(const vec3 &rhs) { return ((*this) = ((*this) * rhs)); }
	inline vec3& operator*=(const T &rhs) { return ((*this) = ((*this) * rhs)); }

	inline vec3& operator/=(const vec3 &rhs) { return ((*this) = ((*this) / rhs)); }
	inline vec3& operator/=(const T &rhs) { return ((*this) = ((*this) / rhs)); }

	inline vec3& operator%=(const vec3 &rhs) { return ((*this) = ((*this) % rhs)); }
	inline vec3& operator%=(const T &rhs) { return ((*this) = ((*this) % rhs)); }

	vec3& operator=(const T &rhs)
	{
		(*this) = vec3(rhs);

		return (*this);
	}

	template<typename T2>
	vec3& operator=(const T2 &rhs)
	{
		(*this) = vec3(rhs);

		return (*this);
	}

	vec3& operator=(const vec3 &rhs)
	{
		this->x = rhs.x;
		this->y = rhs.y;
		this->z = rhs.z;

		return (*this);
	}

	template<typename T2>
	vec3& operator=(const vec3_t<T2> &rhs)
	{
		this->x = T(rhs.x);
		this->y = T(rhs.y);
		this->z = T(rhs.z);

		return (*this);
	}

#pragma endregion

#pragma region Logical Operators
#pragma endregion
#pragma region Comparison Operators

	bool operator==(const vec3 &rhs) const;

	friend inline bool operator==(const vec3 &lhs, const T &rhs) { return (lhs == vec3(rhs)); }
	friend inline bool operator==(const T &lhs, const vec3 &rhs) { return (vec3(lhs) == rhs); }

	friend inline bool operator!=(const vec3 &lhs, const vec3 &rhs) { return !(lhs == rhs); }
	friend inline bool operator!=(const vec3 &lhs, const T &rhs) { return (lhs != vec3(rhs)); }
	friend inline bool operator!=(const T &lhs, const vec3 &rhs) { return (vec3(lhs) != rhs); }

	friend inline bool operator>(const vec3 &lhs, const vec3 &rhs) { return ((lhs.x > rhs.x) && (lhs.y > rhs.y) && (lhs.z > rhs.z)); }
	friend inline bool operator>=(const vec3 &lhs, const vec3 &rhs) { return ((lhs.x >= rhs.x) && (lhs.y >= rhs.y) && (lhs.z >= rhs.z)); }
	friend inline bool operator<(const vec3 &lhs, const vec3 &rhs) { return ((lhs.x < rhs.x) && (lhs.y < rhs.y) && (lhs.z < rhs.z)); }
	friend inline bool operator<=(const vec3 &lhs, const vec3 &rhs) { return ((lhs.x <= rhs.x) && (lhs.y <= rhs.y) && (lhs.z <= rhs.z)); }

	friend inline bool operator>(const vec3 &lhs, const T &rhs) { return (lhs > vec3(rhs)); }
	friend inline bool operator>(const T &lhs, const vec3 &rhs) { return (vec3(lhs) > rhs); }

	friend inline bool operator>=(const vec3 &lhs, const T &rhs) { return (lhs >= vec3(rhs)); }
	friend inline bool operator>=(const T &lhs, const vec3 &rhs) { return (vec3(lhs) >= rhs); }

	friend inline bool operator<(const vec3 &lhs, const T &rhs) { return (lhs < vec3(rhs)); }
	friend inline bool operator<(const T &lhs, const vec3 &rhs) { return (vec3(lhs) < rhs); }

	friend inline bool operator<=(const vec3 &lhs, const T &rhs) { return (lhs <= vec3(rhs)); }
	friend inline bool operator<=(const T &lhs, const vec3 &rhs) { return (vec3(lhs) <= rhs); }

#pragma endregion

#pragma region Cast Operators

	inline explicit operator T*() const
	{
		return reinterpret_cast<T*>(this);
	}

	template<typename T2>
	inline operator vec3_t<T2>() const
	{
		return vec3_t<T2>(
			static_cast<T2>(this->x),
			static_cast<T2>(this->y),
			static_cast<T2>(this->z)
		);
	}
	
#pragma endregion

#pragma region Stream Operators

#ifdef _IOSTREAM_

	friend inline std::ostream& operator<<(std::ostream &stream, const vec3 &rhs)
	{
		return (stream << "vec3(" << rhs.x << ", " << rhs.y << ", " << rhs.z << ")");
	}

	friend inline std::wostream& operator<<(std::wostream &stream, const vec3 &rhs)
	{
		return (stream << L"vec3(" << rhs.x << L", " << rhs.y << L", " << rhs.z << L")");
	}

	friend inline std::istream& operator>>(std::istream &stream, vec3 &rhs)
	{
		rhs = vec3::zero;

		char c;

		_LINALG_IN_FIND_BEGIN();
		stream >> rhs.x;

		_LINALG_IN_FIND_NEXT();
		stream >> rhs.y;

		_LINALG_IN_FIND_NEXT();
		stream >> rhs.z;

		_LINALG_IN_FIND_END();

		return stream;
	}

	friend inline std::wistream& operator>>(std::wistream &stream, vec3 &rhs)
	{
		rhs = vec3::zero;

		wchar_t c;

		_LINALG_IN_FIND_BEGIN();
		stream >> rhs.x;

		_LINALG_IN_FIND_NEXT();
		stream >> rhs.y;

		_LINALG_IN_FIND_NEXT();
		stream >> rhs.z;

		_LINALG_IN_FIND_END();

		return stream;
	}

#endif

#pragma endregion

#pragma endregion

#pragma region

	inline T dot(const vec3 &rhs) const
	{
		return (this->x * rhs.x + this->y * rhs.y + this->z * rhs.z);
	}
	friend inline T dot(const vec3 &lhs, const vec3 &rhs) { return lhs.dot(rhs); }


	vec3 cross(const vec3 &rhs) const
	{
		return vec3(
			((this->y * rhs.z) - (this->z * rhs.y)),
			((this->z * rhs.x) - (this->x * rhs.z)),
			((this->x * rhs.y) - (this->y * rhs.x))
		);
	}
	friend inline vec3 cross(const vec3 &lhs, const vec3 &rhs) { return lhs.cross(rhs); }


	inline T lengthSquared() const { return (this->x * this->x + this->y * this->y + this->z * this->z); }
	friend inline T lengthSquared(const vec3 &lhs) { return lhs.lengthSquared(); }

	inline T length2() const { return (this->x * this->x + this->y * this->y + this->z * this->z); }
	friend inline T length2(const vec3 &lhs) { return lhs.lengthSquared(); }

	inline T length() const { return sqrt(lengthSquared()); }
	friend inline T length(const vec3 &lhs) { return lhs.length(); }


	inline T distanceSquared(const vec3 &rhs) { return ((*this) - rhs).lengthSquared(); }
	friend inline T distanceSquared(const vec3 &lhs, const vec3 &rhs) { return lhs.distanceSquared(rhs); }

	inline T distance(const vec3 &rhs) { return ((*this) - rhs).length(); }
	friend inline T distance(const vec3 &lhs, const vec3 &rhs) { return lhs.distance(rhs); }


	vec3 normalize(const T &to = T(1.0)) const;
	friend inline vec3 normalize(const vec3 &lhs, const T &to = T(1.0)) { return lhs.normalize(to); }


	T angle(const vec3 &rhs) const
	{
		return (dot(rhs) / (length() * rhs.length()));
	}
	friend inline T angle(const vec3 &lhs, const vec3 &rhs) { return lhs.angle(rhs); }


	// Calculates the projection of a onto b
	// 
	// Reference: http://en.wikipedia.org/wiki/Vector_projection#Vector_projection_2
	inline vec3 project(const vec3 &b) const
	{
		const float length = b.length();

		return ((dot(b) / (length * length)) * b);
	}
	friend inline vec3 project(const vec3 &a, const vec3 &b) { return a.project(b); }


	// Calculates the components of a perpendicular to b
	inline vec3 perpendicular(const vec3 &b) const
	{
		const float length = b.length();

		return ((*this) - ((dot(b) / (length * length)) * b));
	}
	friend inline vec3 perpendicular(const vec3 &a, const vec3 &b) { return a.perpendicular(b); }


	// Calculates the reflection vector from entering ray direction a and surface normal b
	inline vec3 reflect(const vec3 &b) const
	{
		return (T(2) * project(b) - (*this));
	}
	friend inline vec3 reflect(const vec3 &a, const vec3 &b) { return a.reflect(b); }


	inline T cosine(const vec3 &b) const
	{
		return normalize().dot(b.normalize());
	}
	friend inline T cosine(const vec3 &a, const vec3 &b) { return a.cosine(b); }

#pragma endregion

#pragma region

	inline bool isNullVector() const;
	friend inline bool isNullVector(const vec3 &v) { return v.isNullVector(); }

	inline bool isUnitVector() const;
	friend inline bool isUnitVector(const vec3 &v) { return v.isUnitVector(); }

	inline bool isNormalized(const T &to = T(1.0)) const;
	friend inline bool isNormalized(const vec3 &v) { return v.isNormalized(); }


	inline bool isOrthogonalTo(const vec3 &rhs) const;
	friend inline bool isOrthogonalTo(const vec3 &lhs, const vec3 &rhs) { return (lhs.isOrthogonalTo(rhs)); }

	inline bool isPerpendicularTo(const vec3 &rhs) const;
	friend inline bool isPerpendicularTo(const vec3 &lhs, const vec3 &rhs) { return (lhs.isPerpendicularTo(rhs)); }

	inline bool isParallelTo(const vec3 &rhs) const;
	friend inline bool isParallelTo(const vec3 &lhs, const vec3 &rhs) { return (lhs.isParallelTo(rhs)); }

#pragma endregion

#pragma region

	inline vec3 abs() const { return vec3(abs(this->x), abs(this->y), abs(this->z)); }
	friend inline vec3 abs(const vec3 &v) { return v.abs(); }


	inline vec3 max(const vec3 &rhs) const
	{
		return vec3(
			((this->x > rhs.x) ? this->x : rhs.x),
			((this->y > rhs.y) ? this->y : rhs.y),
			((this->z > rhs.z) ? this->z : rhs.z)
		);
	}
	friend inline vec3 max(const vec3 &lhs, const vec3 &rhs) { return lhs.max(rhs); }

	inline vec3 min(const vec3 &rhs) const
	{
		return vec3(
			((this->x < rhs.x) ? this->x : rhs.x),
			((this->y < rhs.y) ? this->y : rhs.y),
			((this->z < rhs.z) ? this->z : rhs.z)
		);
	}
	friend inline vec3 min(const vec3 &lhs, const vec3 &rhs) { return lhs.min(rhs); }

	inline vec3 clamp(const vec3 &min, const vec3 &max) const
	{
		return min(max).max(min);
	}
	friend inline vec3 clamp(const vec3 &v, const vec3 &min, const vec3 &max) { return v.clamp(min, max); }


	inline vec3 lerp(const vec3 &to, const T &t) const { return ((*this) + t * (to - (*this))); }
	friend inline vec3 lerp(const vec3 &from, const vec3 &to, const T &t) { return from.lerp(to, t); }

	inline vec3 lerp(const vec3 &to, const vec3 &t) const { return ((*this) + t * (to - (*this))); }
	friend inline vec3 lerp(const vec3 &from, const vec3 &to, const vec3 &t) { return from.lerp(to, t); }


	// Reference: https://en.wikipedia.org/wiki/Slerp
	// Referemce: https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
	inline vec3 slerp(const vec3 &to, const T &t) const
	{
		const T dot = normalize().dot(to.normalize());
		const T theta = acos(dot);
		const T s = sin(theta);

		return (sin((T(1) - t) * theta) / s * (*this) + sin(t * theta) / s * to);
	}
	friend inline vec3 slerp(const vec3 &from, const vec3 &to, const T &t) { return from.slerp(to, t); }

	// Reference: https://en.wikipedia.org/wiki/Slerp
	// Referemce: https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
	inline vec3 slerp(const vec3 &to, const vec3 &t) const
	{
		const T dot = normalize().dot(to.normalize());
		const T theta = acos(dot);
		const T s = sin(theta);

		return (sin((T(1) - t) * theta) / s * (*this) + sin(t * theta) / s * to);
	}
	friend inline vec3 slerp(const vec3 &from, const vec3 &to, const vec3 &t) { return from.slerp(to, t); }


	// Reference: https://en.wikipedia.org/wiki/Sign_function
	inline const vec3 signum() const
	{
		return vec3(
			((x < T(0)) ? T(-1) : ((x > T(0)) ? T(1) : T(0))),
			((y < T(0)) ? T(-1) : ((y > T(0)) ? T(1) : T(0))),
			((z < T(0)) ? T(-1) : ((z > T(0)) ? T(1) : T(0)))
		);
	}
	friend inline vec3 signum(const vec3 &v) { return v.signum(); }

#pragma endregion

#pragma region Swizzling

#define LINALG_SWIZZLE_INDEX(index, c) \
		if ((c == 'X') || (c == 'X') || (c == 'r') || (c == 'R') || (c == 's') || (c == 'S')) index = 0; \
		else if ((c == 'y') || (c == 'Y') || (c == 'g') || (c == 'G') || (c == 't') || (c == 'T')) index = 1; \
		else if ((c == 'z') || (c == 'Z') || (c == 'b') || (c == 'B') || (c == 'p') || (c == 'P')) index = 2; \
		else index = 0;

	vec2 swizzle(const char x, const char y) const
	{
		int x_index = 0, y_index = 1;

		LINALG_SWIZZLE_INDEX(x_index, x);
		LINALG_SWIZZLE_INDEX(y_index, y);

		return vec2((*this)[x_index], (*this)[y_index]);
	}

	vec3 swizzle(const char x, const char y, const char z) const
	{
		int x_index = 0, y_index = 1, z_index = 2;

		LINALG_SWIZZLE_INDEX(x_index, x);
		LINALG_SWIZZLE_INDEX(y_index, y);
		LINALG_SWIZZLE_INDEX(z_index, z);

		return vec3((*this)[x_index], (*this)[y_index], (*this)[z_index]);
	}

#undef LINALG_SWIZZLE_INDEX

#pragma endregion


	inline void swap(vec3 &other)
	{
		const vec3 tmp(*this);
		(*this) = other;
		other = tmp;
	}
	friend inline void swap(vec3 &a, vec3 &b) { a.swap(b); }
};


template<typename T>
class vec4_t
{
private:

	typedef vec2_t<T> vec2;
	typedef vec3_t<T> vec3;
	typedef vec4_t<T> vec4;


public:

	static const vec4_t<T> zero;
	static const vec4_t<T> one;


public:

	T x, y, z, w;


public:

	vec4_t() : x(T(0)), y(T(0)), z(T(0)), w(T(0)) {}

	vec4_t(const vec4_t<T> &v) : x(v.x), y(v.y), z(v.z), w(v.w) {}
	template<typename T2> vec4_t(const vec4_t<T2> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(T(v.w)) {}

	template<typename T2> vec4_t(const T2 &xyzw) : x(T(xyzw)), y(T(xyzw)), z(T(xyzw)), w(T(xyzw)) {}
	template<typename T2> vec4_t(const T2 &x, const T2 &y, const T2 &z, const T2 &w) : x(T(x)), y(T(y)), z(T(z)), w(T(w)) {}

	template<typename T2> vec4_t(const T2 *xyzw) : x(T(xyzw[0])), y(T(xyzw[1])), z(T(xyzw[2])), w(T(xyzw[3])) {}

	template<typename T2, typename T3> vec4_t(const vec2_t<T2> &xy, const T3 &z = T3(0), const T3 &w = T3(0)) : x(T(xy.x)), y(T(xy.y)), z(T(z)), w(T(w)) {}
	template<typename T2, typename T3> vec4_t(const T2 &x, const T2 &y, const vec2_t<T3> &zw) : x(T(x)), y(T(y)), z(T(zw.x)), w(T(zw.y)) {}
	template<typename T2, typename T3> vec4_t(const T2 &x, const vec2_t<T3> &yz, const T2 &w = T3(0)) : x(T(x)), y(T(yz.x)), z(T(yz.y)), w(T(w)) {}
	template<typename T2, typename T3> vec4_t(const vec2_t<T2> &xy, const vec2_t<T3> &zw) : x(T(xy.x)), y(T(xy.y)), z(T(zw.x)), w(T(zw.y)) {}

	template<typename T2, typename T3> vec4_t(const vec3_t<T2> &xyz, const T3 &w = T3(0)) : x(T(xyz.x)), y(T(xyz.y)), z(T(xyz.z)), w(T(w)) {}
	template<typename T2, typename T3> vec4_t(const T2 &x, const vec3_t<T3> &yzw) : x(T(x)), y(T(yzw.x)), z(T(yzw.y)), w(T(yzw.z)) {}

	~vec4_t() {}


#pragma region Operator Overloading

#pragma region Member Access Operators

	inline T& operator[](const int index) { return (reinterpret_cast<T*>(this))[index]; }
	inline T operator[](const int index) const { return ((T*) this)[index]; }

#pragma endregion

#pragma region Arithmetic Operators

	vec4 operator+() const { return vec4(+this->x, +this->y, +this->z, +this->w); }
	vec4 operator-() const { return vec4(-this->x, -this->y, -this->z, -this->w); }

	friend vec4 operator+(const vec4 &lhs, const vec4 &rhs) { return vec4((lhs.x + rhs.x), (lhs.y + rhs.y), (lhs.z + rhs.z), (lhs.w + rhs.w)); }
	friend vec4 operator-(const vec4 &lhs, const vec4 &rhs) { return vec4((lhs.x - rhs.x), (lhs.y - rhs.y), (lhs.z - rhs.z), (lhs.w - rhs.w)); }
	friend vec4 operator*(const vec4 &lhs, const vec4 &rhs) { return vec4((lhs.x * rhs.x), (lhs.y * rhs.y), (lhs.z * rhs.z), (lhs.w * rhs.w)); }
	friend vec4 operator/(const vec4 &lhs, const vec4 &rhs) { return vec4((lhs.x / rhs.x), (lhs.y / rhs.y), (lhs.z / rhs.z), (lhs.w / rhs.w)); }
	friend vec4 operator%(const vec4 &lhs, const vec4 &rhs) { return vec4((lhs.x % rhs.x), (lhs.y % rhs.y), (lhs.z % rhs.z), (lhs.w % rhs.w)); }

	friend inline vec4 operator+(const vec4 &lhs, const T &rhs) { return (lhs + vec4(rhs)); }
	friend inline vec4 operator+(const T &lhs, const vec4 &rhs) { return (vec4(lhs) + rhs); }

	friend inline vec4 operator-(const vec4 &lhs, const T &rhs) { return (lhs - vec4(rhs)); }
	friend inline vec4 operator-(const T &lhs, const vec4 &rhs) { return (vec4(lhs) - rhs); }

	friend inline vec4 operator*(const vec4 &lhs, const T &rhs) { return (lhs * vec4(rhs)); }
	friend inline vec4 operator*(const T &lhs, const vec4 &rhs) { return (vec4(lhs) * rhs); }

	friend inline vec4 operator/(const vec4 &lhs, const T &rhs) { return (lhs / vec4(rhs)); }
	friend inline vec4 operator/(const T &lhs, const vec4 &rhs) { return (vec4(lhs) / rhs); }

	friend inline vec4 operator%(const vec4 &lhs, const T &rhs) { return (lhs % vec4(rhs)); }
	friend inline vec4 operator%(const T &lhs, const vec4 &rhs) { return (vec4(lhs) % rhs); }

#pragma endregion
#pragma region Increment & Decrement Operators

	vec4& operator++() // Prefix
	{
		++this->x;
		++this->y;
		++this->z;
		++this->w;

		return (*this);
	}
	inline vec4 operator++(int) { return ++(*this); } // Postfix

	vec4& operator--() // Prefix
	{
		--this->x;
		--this->y;
		--this->z;
		--this->w;

		return (*this);
	}
	inline vec4 operator--(int) { return --(*this); } // Postfix

#pragma endregion
#pragma region Assignment Operators

	inline vec4& operator+=(const vec4 &rhs) { return ((*this) = ((*this) + rhs)); }
	inline vec4& operator+=(const T &rhs) { return ((*this) = ((*this) + rhs)); }

	inline vec4& operator-=(const vec4 &rhs) { return ((*this) = ((*this) - rhs)); }
	inline vec4& operator-=(const T &rhs) { return ((*this) = ((*this) - rhs)); }

	inline vec4& operator*=(const vec4 &rhs) { return ((*this) = ((*this) * rhs)); }
	inline vec4& operator*=(const T &rhs) { return ((*this) = ((*this) * rhs)); }

	inline vec4& operator/=(const vec4 &rhs) { return ((*this) = ((*this) / rhs)); }
	inline vec4& operator/=(const T &rhs) { return ((*this) = ((*this) / rhs)); }

	inline vec4& operator%=(const vec4 &rhs) { return ((*this) = ((*this) % rhs)); }
	inline vec4& operator%=(const T &rhs) { return ((*this) = ((*this) % rhs)); }

	vec4& operator=(const T &rhs)
	{
		(*this) = vec4(rhs);

		return (*this);
	}

	template<typename T2>
	vec4& operator=(const T2 &rhs)
	{
		(*this) = vec4(rhs);

		return (*this);
	}

	vec4& operator=(const vec4 &rhs)
	{
		this->x = rhs.x;
		this->y = rhs.y;
		this->z = rhs.z;
		this->w = rhs.w;

		return (*this);
	}

	template<typename T2>
	vec4& operator=(const vec4_t<T2> &rhs)
	{
		this->x = T(rhs.x);
		this->y = T(rhs.y);
		this->z = T(rhs.z);
		this->w = T(rhs.w);

		return (*this);
	}

#pragma endregion

#pragma region Logical Operators
#pragma endregion
#pragma region Comparison Operators

	bool operator==(const vec4 &rhs) const;

	friend inline bool operator==(const vec4 &lhs, const T &rhs) { return (lhs == vec4(rhs)); }
	friend inline bool operator==(const T &lhs, const vec4 &rhs) { return (vec4(lhs) == rhs); }

	friend inline bool operator!=(const vec4 &lhs, const vec4 &rhs) { return !(lhs == rhs); }
	friend inline bool operator!=(const vec4 &lhs, const T &rhs) { return (lhs != vec4(rhs)); }
	friend inline bool operator!=(const T &lhs, const vec4 &rhs) { return (vec4(lhs) != rhs); }

	friend inline bool operator>(const vec4 &lhs, const vec4 &rhs) { return ((lhs.x > rhs.x) && (lhs.y > rhs.y) && (lhs.z > rhs.z) && (lhs.w > rhs.w)); }
	friend inline bool operator>=(const vec4 &lhs, const vec4 &rhs) { return ((lhs.x >= rhs.x) && (lhs.y >= rhs.y) && (lhs.z >= rhs.z) && (lhs.w >= rhs.w)); }
	friend inline bool operator<(const vec4 &lhs, const vec4 &rhs) { return ((lhs.x < rhs.x) && (lhs.y < rhs.y) && (lhs.z < rhs.z) && (lhs.w < rhs.w)); }
	friend inline bool operator<=(const vec4 &lhs, const vec4 &rhs) { return ((lhs.x <= rhs.x) && (lhs.y <= rhs.y) && (lhs.z <= rhs.z) && (lhs.w <= rhs.w)); }

	friend inline bool operator>(const vec4 &lhs, const T &rhs) { return (lhs > vec4(rhs)); }
	friend inline bool operator>(const T &lhs, const vec4 &rhs) { return (vec4(lhs) > rhs); }

	friend inline bool operator>=(const vec4 &lhs, const T &rhs) { return (lhs >= vec4(rhs)); }
	friend inline bool operator>=(const T &lhs, const vec4 &rhs) { return (vec4(lhs) >= rhs); }

	friend inline bool operator<(const vec4 &lhs, const T &rhs) { return (lhs < vec4(rhs)); }
	friend inline bool operator<(const T &lhs, const vec4 &rhs) { return (vec4(lhs) < rhs); }

	friend inline bool operator<=(const vec4 &lhs, const T &rhs) { return (lhs <= vec4(rhs)); }
	friend inline bool operator<=(const T &lhs, const vec4 &rhs) { return (vec4(lhs) <= rhs); }

#pragma endregion

#pragma region Cast Operators

	inline explicit operator T*() const
	{
		return reinterpret_cast<T*>(this);
	}

	template<typename T2>
	inline operator vec4_t<T2>() const
	{
		return vec4_t<T2>(
			static_cast<T2>(this->x),
			static_cast<T2>(this->y),
			static_cast<T2>(this->z),
			static_cast<T2>(this->w)
		);
	}
	
#pragma endregion

#pragma region Stream Operators

#ifdef _IOSTREAM_

	friend inline std::ostream& operator<<(std::ostream &stream, const vec4 &rhs)
	{
		return (stream << "vec4(" << rhs.x << ", " << rhs.y << ", " << rhs.z << ", " << rhs.w << ")");
	}

	friend inline std::wostream& operator<<(std::wostream &stream, const vec4 &rhs)
	{
		return (stream << L"vec4(" << rhs.x << L", " << rhs.y << L", " << rhs.z << L", " << rhs.w << L")");
	}

	friend inline std::istream& operator>>(std::istream &stream, vec4 &rhs)
	{
		rhs = vec4::zero;

		char c;

		_LINALG_IN_FIND_BEGIN();
		stream >> rhs.x;

		_LINALG_IN_FIND_NEXT();
		stream >> rhs.y;

		_LINALG_IN_FIND_NEXT();
		stream >> rhs.z;

		_LINALG_IN_FIND_NEXT();
		stream >> rhs.w;

		_LINALG_IN_FIND_END();

		return stream;
	}

	friend inline std::wistream& operator>>(std::wistream &stream, vec4 &rhs)
	{
		rhs = vec4::zero;

		wchar_t c;

		_LINALG_IN_FIND_BEGIN();
		stream >> rhs.x;

		_LINALG_IN_FIND_NEXT();
		stream >> rhs.y;

		_LINALG_IN_FIND_NEXT();
		stream >> rhs.z;

		_LINALG_IN_FIND_NEXT();
		stream >> rhs.w;

		_LINALG_IN_FIND_END();

		return stream;
	}

#endif

#pragma endregion

#pragma endregion

#pragma region

	inline T dot(const vec4 &rhs) const
	{
		return (this->x * rhs.x + this->y * rhs.y + this->z * rhs.z + this->w * rhs.w);
	}
	friend inline T dot(const vec4 &lhs, const vec4 &rhs) { return lhs.dot(rhs); }


	inline T lengthSquared() const { return (this->x * this->x + this->y * this->y + this->z * this->z + this->w * this->w); }
	friend inline T lengthSquared(const vec4 &lhs) { return lhs.lengthSquared(); }

	inline T length2() const { return (this->x * this->x + this->y * this->y + this->z * this->z + this->w * this->w); }
	friend inline T length2(const vec4 &lhs) { return lhs.lengthSquared(); }

	inline T length() const { return sqrt(lengthSquared()); }
	friend inline T length(const vec4 &lhs) { return lhs.length(); }


	inline T distanceSquared(const vec4 &rhs) { return ((*this) - rhs).lengthSquared(); }
	friend inline T distanceSquared(const vec4 &lhs, const vec4 &rhs) { return lhs.distanceSquared(rhs); }

	inline T distance(const vec4 &rhs) { return ((*this) - rhs).length(); }
	friend inline T distance(const vec4 &lhs, const vec4 &rhs) { return lhs.distance(rhs); }


	vec4 normalize(const T &to = T(1.0)) const;
	friend inline vec4 normalize(const vec4 &lhs, const T &to = T(1.0)) { return lhs.normalize(to); }


	T angle(const vec4 &rhs) const
	{
		return (dot(rhs) / (length() * rhs.length()));
	}
	friend inline T angle(const vec4 &lhs, const vec4 &rhs) { return lhs.angle(rhs); }

#pragma endregion

#pragma region

	inline bool isNullVector() const;
	friend inline bool isNullVector(const vec4 &v) { return v.isNullVector(); }

	inline bool isUnitVector() const;
	friend inline bool isUnitVector(const vec4 &v) { return v.isUnitVector(); }

	inline bool isNormalized(const T &to = T(1.0)) const;
	friend inline bool isNormalized(const vec4 &v) { return v.isNormalized(); }


	inline bool isOrthogonalTo(const vec4 &rhs) const;
	friend inline bool isOrthogonalTo(const vec4 &lhs, const vec4 &rhs) { return (lhs.isOrthogonalTo(rhs)); }

	inline bool isPerpendicularTo(const vec4 &rhs) const;
	friend inline bool isPerpendicularTo(const vec4 &lhs, const vec4 &rhs) { return (lhs.isPerpendicularTo(rhs)); }

	inline bool isParallelTo(const vec4 &rhs) const;
	friend inline bool isParallelTo(const vec4 &lhs, const vec4 &rhs) { return (lhs.isParallelTo(rhs)); }

#pragma endregion

#pragma region

	inline vec4 abs() const { return vec4(abs(this->x), abs(this->y), abs(this->z), abs(this->w)); }
	friend inline vec4 abs(const vec4 &v) { return v.abs(); }


	inline vec4 max(const vec4 &rhs) const
	{
		return vec4(
			((this->x > rhs.x) ? this->x : rhs.x),
			((this->y > rhs.y) ? this->y : rhs.y),
			((this->z > rhs.z) ? this->z : rhs.z),
			((this->w > rhs.w) ? this->w : rhs.w)
		);
	}
	friend inline vec4 max(const vec4 &lhs, const vec4 &rhs) { return lhs.max(rhs); }

	inline vec4 min(const vec4 &rhs) const
	{
		return vec4(
			((this->x < rhs.x) ? this->x : rhs.x),
			((this->y < rhs.y) ? this->y : rhs.y),
			((this->z < rhs.z) ? this->z : rhs.z),
			((this->w < rhs.w) ? this->w : rhs.w)
		);
	}
	friend inline vec4 min(const vec4 &lhs, const vec4 &rhs) { return lhs.min(rhs); }

	inline vec4 clamp(const vec4 &min, const vec4 &max) const
	{
		return min(max).max(min);
	}
	friend inline vec4 clamp(const vec4 &v, const vec4 &min, const vec4 &max) { return v.clamp(min, max); }


	inline vec4 lerp(const vec4 &to, const T &t) const { return ((*this) + t * (to - (*this))); }
	friend inline vec4 lerp(const vec4 &from, const vec4 &to, const T &t) { return from.lerp(to, t); }

	inline vec4 lerp(const vec4 &to, const vec4 &t) const { return ((*this) + t * (to - (*this))); }
	friend inline vec4 lerp(const vec4 &from, const vec4 &to, const vec4 &t) { return from.lerp(to, t); }


	// Reference: https://en.wikipedia.org/wiki/Slerp
	// Referemce: https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
	vec4 slerp(const vec4 &to, const T &t) const
	{
		const T dot = normalize().dot(to.normalize());
		const T theta = acos(dot);
		const T s = sin(theta);

		return (sin((T(1) - t) * theta) / s * (*this) + sin(t * theta) / s * to);
	}
	friend inline vec4 slerp(const vec4 &from, const vec4 &to, const T &t) { return from.slerp(to, t); }

	// Reference: https://en.wikipedia.org/wiki/Slerp
	// Referemce: https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
	vec4 slerp(const vec4 &to, const vec4 &t) const
	{
		const T dot = normalize().dot(to.normalize());
		const T theta = acos(dot);
		const T s = sin(theta);

		return (sin((T(1) - t) * theta) / s * (*this) + sin(t * theta) / s * to);
	}
	friend inline vec4 slerp(const vec4 &from, const vec4 &to, const vec4 &t) { return from.slerp(to, t); }


	// Reference: https://en.wikipedia.org/wiki/Sign_function
	inline const vec4 signum() const
	{
		return vec4(
			((this->x < T(0)) ? T(-1) : ((this->x > T(0)) ? T(1) : T(0))),
			((this->y < T(0)) ? T(-1) : ((this->y > T(0)) ? T(1) : T(0))),
			((this->z < T(0)) ? T(-1) : ((this->z > T(0)) ? T(1) : T(0))),
			((this->w < T(0)) ? T(-1) : ((this->w > T(0)) ? T(1) : T(0)))
		);
	}
	friend inline vec4 signum(const vec4 &v) { return v.signum(); }

#pragma endregion

#pragma region Swizzling

#define LINALG_SWIZZLE_INDEX(index, c) \
		if ((c == 'X') || (c == 'X') || (c == 'r') || (c == 'R') || (c == 's') || (c == 'S')) index = 0; \
		else if ((c == 'y') || (c == 'Y') || (c == 'g') || (c == 'G') || (c == 't') || (c == 'T')) index = 1; \
		else if ((c == 'z') || (c == 'Z') || (c == 'b') || (c == 'B') || (c == 'p') || (c == 'P')) index = 2; \
		else if ((c == 'w') || (c == 'W') || (c == 'a') || (c == 'A') || (c == 'q') || (c == 'Q')) index = 3; \
		else index = 0;

	vec2 swizzle(const char x, const char y) const
	{
		int x_index = 0, y_index = 1;

		LINALG_SWIZZLE_INDEX(x_index, x);
		LINALG_SWIZZLE_INDEX(y_index, y);

		return vec2((*this)[x_index], (*this)[y_index]);
	}

	vec3 swizzle(const char x, const char y, const char z) const
	{
		int x_index = 0, y_index = 1, z_index = 2;

		LINALG_SWIZZLE_INDEX(x_index, x);
		LINALG_SWIZZLE_INDEX(y_index, y);
		LINALG_SWIZZLE_INDEX(z_index, z);

		return vec3((*this)[x_index], (*this)[y_index], (*this)[z_index]);
	}

	vec4 swizzle(const char x, const char y, const char z, const char w) const
	{
		int x_index = 0, y_index = 1, z_index = 2, w_index = 3;

		LINALG_SWIZZLE_INDEX(x_index, x);
		LINALG_SWIZZLE_INDEX(y_index, y);
		LINALG_SWIZZLE_INDEX(z_index, z);
		LINALG_SWIZZLE_INDEX(w_index, w);

		return vec4((*this)[x_index], (*this)[y_index], (*this)[z_index], (*this)[w_index]);
	}

#undef LINALG_SWIZZLE_INDEX

#pragma endregion


	inline void swap(vec4 &other)
	{
		const vec4 tmp(*this);
		(*this) = other;
		other = tmp;
	}
	friend inline void swap(vec4 &a, vec4 &b) { a.swap(b); }
};


template<typename T>
class mat2_t
{
private:

	typedef vec2_t<T> vec2;

	typedef mat2_t<T> mat2;
	typedef mat3_t<T> mat3;
	typedef mat4_t<T> mat4;


public:

	static const mat2_t<T> zero;
	static const mat2_t<T> identity;


public:

	vec2 columns[2];


public:

	mat2_t(const T mainDiagonalValue = T(1))
	{
		this->columns[0] = vec2(mainDiagonalValue, T(0));
		this->columns[1] = vec2(T(0), mainDiagonalValue);
	}

	mat2_t(
		const vec2 &column1, // first column
		const vec2 &column2) // second column
	{
		this->columns[0] = column1;
		this->columns[1] = column2;
	}

	mat2_t(const vec2 columns[2])
	{
		this->columns[0] = columns[0];
		this->columns[1] = columns[1];
	}

	mat2_t(const T values[2 * 2])
	{
		for (int i = 0; i < 2 * 2; i++)
			(reinterpret_cast<T*>(this))[i] = values[i];
	}

	mat2_t(
		const T a, const T b, // first column
		const T c, const T d) // second column
	{
		(*this) = mat2(
			vec2(a, b),
			vec2(c, d)
		);
	}

	template<typename T2> mat2_t(const mat3_t<T2> &m);
	template<typename T2> mat2_t(const mat4_t<T2> &m);

	~mat2_t() {}


#pragma region Operator Overloading

#pragma region Member Access Operators

	inline vec2& operator[](const int index) { return (reinterpret_cast<vec2*>(this))[index]; }
	inline vec2 operator[](const int index) const { return ((vec2*) this)[index]; }

#pragma endregion

#pragma region Arithmetic Operators

	mat2 operator+(const mat2 &rhs) const
	{
		mat2 result;

		for (int i = 0; i < 2; i++)
			result[i] = (*this)[i] + rhs[i];

		return result;
	}

	mat2 operator-(const mat2 &rhs) const
	{
		mat2 result;

		for (int i = 0; i < 2; i++)
			result[i] = (*this)[i] - rhs[i];

		return result;
	}

	mat2 operator*(const mat2 &rhs) const
	{
		mat2 result;

		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				result.value(i, j, row(i).dot(rhs.col(j)));

		return result;
	}

	vec2 operator*(const vec2 &rhs) const
	{
		return vec2(
			(rhs.x * (*this)[0].x) + (rhs.y * (*this)[1].x),
			(rhs.x * (*this)[0].y) + (rhs.y * (*this)[1].y)
		);
	}

	friend vec2 operator*(const vec2 &lhs, const mat2 &rhs)
	{
		return vec2(
			lhs.dot(rhs[0]),
			lhs.dot(rhs[1])
		);
	}

	mat2 operator*(const T &rhs) const
	{
		mat2 result;

		for (int i = 0; i < 2; i++)
			result[i] = (*this)[i] * rhs;

		return result;
	}
	friend inline mat2 operator*(const T &lhs, const mat2 &rhs) { return (rhs * lhs); }

	mat2 operator/(const T &rhs) const
	{
		mat2 result;

		for (int i = 0; i < 2; i++)
			result[i] = (*this)[i] / rhs;

		return result;
	}

#pragma endregion
#pragma region Assignment Operators

	mat2& operator+=(const mat2 &rhs) { return ((*this) = (*this) + rhs); }
	mat2& operator-=(const mat2 &rhs) { return ((*this) = (*this) - rhs); }
	mat2& operator*=(const mat2 &rhs) { return ((*this) = (*this) * rhs); }
	mat2& operator*=(const T rhs) { return ((*this) = (*this) * rhs); }
	mat2& operator/=(const T rhs) { return ((*this) = (*this) / rhs); }

	mat2& operator=(const T &rhs)
	{
		(*this) = mat2(rhs);

		return (*this);
	}

	template<typename T2>
	mat2& operator=(const T2 &rhs)
	{
		(*this) = mat2(rhs);

		return (*this);
	}

	mat2& operator=(const mat2 &rhs)
	{
		for (int i = 0; i < 2; i++)
			this->columns[i] = rhs[i];

		return (*this);
	}

	template<typename T2>
	mat2& operator=(const mat2_t<T2> &rhs)
	{
		for (int i = 0; i < 2; i++)
			this->columns[i] = vec2(rhs[i]);

		return (*this);
	}

#pragma endregion

#pragma region Comparison Operators
	
	bool operator==(const mat2 &rhs) const;

	inline bool operator!=(const mat2 &rhs) const { return !((*this) == rhs); }

#pragma endregion

#pragma region Cast Operators

	inline explicit operator T*() const
	{
		return reinterpret_cast<T*>(this);
	}

	template<typename T2>
	inline operator mat2_t<T2>() const
	{
		return mat2_t<T2>(
			static_cast<vec2_t<T2>>(this->columns[0]),
			static_cast<vec2_t<T2>>(this->columns[1])
		);
	}

#pragma endregion

#pragma region Stream Operators

#ifdef _IOSTREAM_

	friend inline std::ostream& operator<<(std::ostream &stream, const mat2 &rhs)
	{
		return (stream << "mat2(" << rhs[0] << "," << std::endl
					   << "     " << rhs[1] << ")");
	}

	friend inline std::wostream& operator<<(std::wostream &stream, const mat2 &rhs)
	{
		return (stream << L"mat2(" << rhs[0] << L"," << std::endl
					   << L"     " << rhs[1] << L")");
	}

	friend inline std::istream& operator>>(std::istream &stream, mat2 &rhs)
	{
		rhs = mat2::identity;

		char c;

		_LINALG_IN_FIND_BEGIN();
		stream >> rhs[0];

		_LINALG_IN_FIND_NEXT();
		stream >> rhs[1];

		_LINALG_IN_FIND_END();

		return stream;
	}

	friend inline std::wistream& operator>>(std::wistream &stream, mat2 &rhs)
	{
		rhs = mat2::identity;

		wchar_t c;

		_LINALG_IN_FIND_BEGIN();
		stream >> rhs[0];

		_LINALG_IN_FIND_NEXT();
		stream >> rhs[1];

		_LINALG_IN_FIND_END();

		return stream;
	}

#endif

#pragma endregion

#pragma endregion


	T determinant() const
	{
		return ((*this)[0][0] * (*this)[1][1]) - ((*this)[1][0] * (*this)[0][1]);
	}
	friend inline T determinant(const mat2 &m) { return mat2(m).determinant(); }


	mat2& inverse()
	{
		const T d = T(1) / determinant();

		return ((*this) = (d * mat2(
			(*this)[1][1], -(*this)[0][1],
			-(*this)[1][0], (*this)[0][0]
		)));
	}
	friend inline mat2 inverse(const mat2 &m) { return mat2(m).inverse(); }


	mat2& transpose()
	{
		return ((*this) = mat2(
			vec2((*this)[0].x, (*this)[1].x),
			vec2((*this)[0].y, (*this)[1].y)
		));
	}
	friend inline mat2 transpose(const mat2 &m) { return mat2(m).transpose(); }


	inline vec2 col(const int index) const
	{
		return (*this)[index];
	}

	inline vec2 row(const int index) const
	{
		return vec2(
			(*this)[0][index],
			(*this)[1][index]
		);
	}


	inline T value(const int row, const int column) const
	{
		return (*this)[column][row];
	}

	inline void value(const int row, const int column, const T &value)
	{
		(*this)[column][row] = value;
	}


	inline mat2& setZero()
	{
		return ((*this) = mat2::zero);
	}

	inline mat2& setIdentity()
	{
		return ((*this) = mat2::identity);
	}


	inline void swap(mat2 &other)
	{
		const mat2 tmp(*this);
		(*this) = other;
		other = tmp;
	}
	friend inline void swap(mat2 &a, mat2 &b) { a.swap(b); }
};


template<typename T>
class mat3_t
{
private:

	typedef vec2_t<T> vec2;
	typedef vec3_t<T> vec3;

	typedef mat2_t<T> mat2;
	typedef mat3_t<T> mat3;
	typedef mat4_t<T> mat4;


public:

	static const mat3_t<T> zero;
	static const mat3_t<T> identity;


public:

	vec3 columns[3];


public:

	mat3_t(const T mainDiagonalValue = T(1))
	{
		this->columns[0] = vec3(mainDiagonalValue, T(0), T(0));
		this->columns[1] = vec3(T(0), mainDiagonalValue, T(0));
		this->columns[2] = vec3(T(0), T(0), mainDiagonalValue);
	}

	mat3_t(
		const vec3 &column1, // first column
		const vec3 &column2, // second column
		const vec3 &column3) // third column
	{
		this->columns[0] = column1;
		this->columns[1] = column2;
		this->columns[2] = column3;
	}

	mat3_t(const vec3 columns[3])
	{
		this->columns[0] = columns[0];
		this->columns[1] = columns[1];
		this->columns[2] = columns[2];
	}

	mat3_t(const T values[3 * 3])
	{
		for (int i = 0; i < 3 * 3; i++)
			(reinterpret_cast<T*>(this))[i] = values[i];
	}

	mat3_t(
		const T a, const T b, const T c, // first column
		const T d, const T e, const T f, // second column
		const T g, const T h, const T i) // third column
	{
		(*this) = mat3(
			vec3(a, b, c),
			vec3(d, e, f),
			vec3(g, h, i)
		);
	}

	template<typename T2> mat3_t(const mat4_t<T2> &m);

	~mat3_t() {}


#pragma region Operator Overloading

#pragma region Member Access Operators

	inline vec3& operator[](const int index) { return (reinterpret_cast<vec3*>(this))[index]; }
	inline vec3 operator[](const int index) const { return ((vec3*) this)[index]; }

#pragma endregion

#pragma region Function Call Operator

	T operator()(const int row, const int column) const
	{
		return (*this)[column][row];
	}

	T& operator()(const int row, const int column)
	{
		return (*this)[column][row];
	}

#pragma endregion

#pragma region Arithmetic Operators

	mat3 operator+(const mat3 &rhs) const
	{
		mat3 result;

		for (int i = 0; i < 3; i++)
			result[i] = (*this)[i] + rhs[i];

		return result;
	}

	mat3 operator-(const mat3 &rhs) const
	{
		mat3 result;

		for (int i = 0; i < 3; i++)
			result[i] = (*this)[i] - rhs[i];

		return result;
	}

	mat3 operator*(const mat3 &rhs) const
	{
		mat3 result;

		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				result.value(i, j, row(i).dot(rhs.col(j)));

		return result;
	}

	vec3 operator*(const vec3 &rhs) const
	{
		return vec3(
			(rhs.x * (*this)[0].x) + (rhs.y * (*this)[1].x) + (rhs.z * (*this)[2].x),
			(rhs.x * (*this)[0].y) + (rhs.y * (*this)[1].y) + (rhs.z * (*this)[2].y),
			(rhs.x * (*this)[0].z) + (rhs.y * (*this)[1].z) + (rhs.z * (*this)[2].z)
		);
	}

	friend vec3 operator*(const vec3 &lhs, const mat3 &rhs)
	{
		return vec3(
			lhs.dot(rhs[0]),
			lhs.dot(rhs[1]),
			lhs.dot(rhs[2])
		);
	}

	mat3 operator*(const T &rhs) const
	{
		mat3 result;

		for (int i = 0; i < 3; i++)
			result[i] = (*this)[i] * rhs;

		return result;
	}
	friend inline mat3 operator*(const T &lhs, const mat3 &rhs) { return (rhs * lhs); }

	mat3 operator/(const T &rhs) const
	{
		mat3 result;

		for (int i = 0; i < 3; i++)
			result[i] = (*this)[i] / rhs;

		return result;
	}

#pragma endregion
#pragma region Assignment Operators

	mat3& operator+=(const mat3 &rhs) { return ((*this) = (*this) + rhs); }
	mat3& operator-=(const mat3 &rhs) { return ((*this) = (*this) - rhs); }
	mat3& operator*=(const mat3 &rhs) { return ((*this) = (*this) * rhs); }
	mat3& operator*=(const T rhs) { return ((*this) = (*this) * rhs); }
	mat3& operator/=(const T rhs) { return ((*this) = (*this) / rhs); }

	mat3& operator=(const T &rhs)
	{
		(*this) = mat3(rhs);

		return (*this);
	}

	template<typename T2>
	mat3& operator=(const T2 &rhs)
	{
		(*this) = mat3(rhs);

		return (*this);
	}

	mat3& operator=(const mat3 &rhs)
	{
		for (int i = 0; i < 3; i++)
			this->columns[i] = rhs[i];

		return (*this);
	}

	template<typename T2>
	mat3& operator=(const mat3_t<T2> &rhs)
	{
		for (int i = 0; i < 3; i++)
			this->columns[i] = vec3(rhs[i]);

		return (*this);
	}

#pragma endregion

#pragma region Comparison Operators

	bool operator==(const mat3 &rhs) const;

	inline bool operator!=(const mat3 &rhs) const { return !((*this) == rhs); }

#pragma endregion

#pragma region Cast Operators

	inline explicit operator T*() const
	{
		return reinterpret_cast<T*>(this);
	}

	template<typename T2>
	inline operator mat3_t<T2>() const
	{
		return mat3_t<T2>(
			static_cast<vec3_t<T2>>(this->columns[0]),
			static_cast<vec3_t<T2>>(this->columns[1]),
			static_cast<vec3_t<T2>>(this->columns[2])
		);
	}

#pragma endregion

#pragma region Stream Operators

#ifdef _IOSTREAM_

	friend inline std::ostream& operator<<(std::ostream &stream, const mat3 &rhs)
	{
		return (stream << "mat3(" << rhs[0] << "," << std::endl
					   << "     " << rhs[1] << "," << std::endl
					   << "     " << rhs[2] << ")");
	}

	friend inline std::wostream& operator<<(std::wostream &stream, const mat3 &rhs)
	{
		return (stream << L"mat3(" << rhs[0] << L"," << std::endl
					   << L"     " << rhs[1] << L"," << std::endl
					   << L"     " << rhs[2] << L")");
	}

	friend inline std::istream& operator>>(std::istream &stream, mat3 &rhs)
	{
		rhs = mat3::identity;

		char c;

		_LINALG_IN_FIND_BEGIN();
		stream >> rhs[0];

		_LINALG_IN_FIND_NEXT();
		stream >> rhs[1];

		_LINALG_IN_FIND_NEXT();
		stream >> rhs[2];

		_LINALG_IN_FIND_END();

		return stream;
	}

	friend inline std::wistream& operator>>(std::wistream &stream, mat3 &rhs)
	{
		rhs = mat3::identity;

		wchar_t c;

		_LINALG_IN_FIND_BEGIN();
		stream >> rhs[0];

		_LINALG_IN_FIND_NEXT();
		stream >> rhs[1];

		_LINALG_IN_FIND_NEXT();
		stream >> rhs[2];

		_LINALG_IN_FIND_END();

		return stream;
	}

#endif

#pragma endregion

#pragma endregion


	T determinant() const
	{
		return (*this)[0][0] * ((*this)[1][1] * (*this)[2][2]) - ((*this)[2][1] * (*this)[1][2])
			- (*this)[1][0] * ((*this)[0][1] * (*this)[2][2]) - ((*this)[2][1] * (*this)[0][2])
			+ (*this)[2][0] * ((*this)[0][1] * (*this)[1][2]) - ((*this)[1][1] * (*this)[0][2]);
	}
	friend inline T determinant(const mat3 &m) { return mat3(m).determinant(); }


	mat3& inverse()
	{
		const T d = T(1) / determinant();

		return ((*this) = (d * mat3(
			((*this)[2][2] * (*this)[1][1]) - ((*this)[1][2] * (*this)[2][1]),
			((*this)[1][2] * (*this)[2][0]) - ((*this)[2][2] * (*this)[1][0]),
			((*this)[2][1] * (*this)[1][0]) - ((*this)[1][1] * (*this)[2][0]),

			((*this)[0][2] * (*this)[2][1]) - ((*this)[2][2] * (*this)[0][1]),
			((*this)[2][2] * (*this)[0][0]) - ((*this)[0][2] * (*this)[2][0]),
			((*this)[0][1] * (*this)[2][0]) - ((*this)[2][1] * (*this)[0][0]),

			((*this)[1][2] * (*this)[0][1]) - ((*this)[0][2] * (*this)[1][1]),
			((*this)[0][2] * (*this)[1][0]) - ((*this)[1][2] * (*this)[0][0]),
			((*this)[1][1] * (*this)[0][0]) - ((*this)[0][1] * (*this)[1][0])
		)));
	}
	friend inline mat3 inverse(const mat3 &m) { return mat3(m).inverse(); }


	mat3& transpose()
	{
		return ((*this) = mat3(
			vec3((*this)[0].x, (*this)[1].x, (*this)[2].x),
			vec3((*this)[0].y, (*this)[1].y, (*this)[2].y),
			vec3((*this)[0].z, (*this)[1].z, (*this)[2].z)
		));
	}
	friend inline mat3 transpose(const mat3 &m) { return mat3(m).transpose(); }


	inline vec3 col(const int index) const
	{
		return (*this)[index];
	}

	inline vec3 row(const int index) const
	{
		return vec3(
			(*this)[0][index],
			(*this)[1][index],
			(*this)[2][index]
		);
	}


	inline T value(const int row, const int column) const
	{
		return (*this)[column][row];
	}

	inline void value(const int row, const int column, const T &value)
	{
		(*this)[column][row] = value;
	}


	inline mat3& setZero()
	{
		return ((*this) = mat3::zero);
	}

	inline mat3& setIdentity()
	{
		return ((*this) = mat3::identity);
	}


	mat3& rotate(const T radians, const vec3 &axis)
	{
		vec3 axisNormalized = axis;

		if (!axisNormalized.isUnitVector())
			axisNormalized = normalize(axisNormalized);

		const T s = sin(radians);
		const T c = cos(radians);
		const T oc = T(1) - c;

		const T x = axisNormalized.x;
		const T y = axisNormalized.y;
		const T z = axisNormalized.z;

		return ((*this) *= mat3(
			(x * x * oc + c),
			(x * y * oc - z * s),
			(x * z * oc + y * s),

			(y * x * oc + z * s),
			(y * y * oc + c),
			(y * z * oc - x * s),

			(x * z * oc - y * s),
			(y * z * oc + x * s),
			(z * z * oc + c)
		));
	}
	friend inline mat3 rotate(const mat3 &m, const T radians, const vec3 &axis) { return mat3(m).rotate(radians, axis); }

	inline mat3& rotate(const T radians, const T ax, const T ay, const T az) { return (*this).rotate(radians, vec3(ax, ay, az)); }
	friend inline mat3 rotate(const mat3 &m, const T radians, const T ax, const T ay, const T az) { return mat3(m).rotate(radians, ax, ay, az); }


	inline mat3& rotateDegrees(const T degrees, const vec3 &axis) { return rotate(degrees * T(LINALG_DEG2RAD), axis); }
	friend inline mat3 rotateDegrees(const mat3 &m, const T degrees, const vec3 &axis) { return mat3(m).rotateDegrees(degrees, axis); }

	inline mat3& rotateDegrees(const T degrees, const T ax, const T ay, const T az) { return (*this).rotateDegrees(degrees, vec3(ax, ay, az)); }
	friend inline mat3 rotateDegrees(const mat3 &m, const T degrees, const T ax, const T ay, const T az) { return mat3(m).rotateDegrees(degrees, ax, ay, az); }


	mat3& rotateX(const T radians)
	{
		const T s = sin(radians);
		const T c = cos(radians);

		return ((*this) *= mat3(
			T(1), T(0), T(0),
			T(0), c, -s,
			T(0), s, c
		));
	}
	friend inline mat3 rotateX(const mat3 &m, const T radians) { return mat3(m).rotateX(radians); }

	inline mat3& rotateXDegrees(const T degrees) { return rotateX(degrees * T(LINALG_DEG2RAD)); }
	friend inline mat3 rotateXDegrees(const mat3 &m, const T degrees) { return mat3(m).rotateXDegrees(degrees); }


	mat3& rotateY(const T radians)
	{
		const T s = sin(radians);
		const T c = cos(radians);

		return ((*this) *= mat3(
			c, T(0), s,
			T(0), T(1), T(0),
			-s, T(0), c
		));
	}
	friend inline mat3 rotateY(const mat3 &m, const T radians) { return mat3(m).rotateY(radians); }

	inline mat3& rotateYDegrees(const T degrees) { return rotateY(degrees * T(LINALG_DEG2RAD)); }
	friend inline mat3 rotateYDegrees(const mat3 &m, const T degrees) { return mat3(m).rotateYDegrees(degrees); }


	mat3& rotateZ(const T radians)
	{
		const T s = sin(radians);
		const T c = cos(radians);

		return ((*this) *= mat3(
			c, -s, T(0),
			s, c, T(0),
			T(0), T(0), T(1)
		));
	}
	friend inline mat3 rotateZ(const mat3 &m, const T radians) { return mat3(m).rotateZ(radians); }

	inline mat3& rotateZDegrees(const T degrees) { return rotateZ(degrees * T(LINALG_DEG2RAD)); }
	friend inline mat3 rotateZDegrees(const mat3 &m, const T degrees) { return mat3(m).rotateZDegrees(degrees); }


	inline mat3& setScaling(const vec3 &scaling)
	{
		return ((*this) = mat3::scaling(scaling));
	}

	inline mat3 setScaling(const T sx, const T sy, const T sz = T(0))
	{
		return ((*this) = mat3::scaling(vec3(sx, sy, sz)));
	}

	inline vec3 getScaling() const
	{
		return vec3((*this)(0, 0), (*this)(1, 1), (*this)(2, 2));
	}


	// TODO: Gives incorrect results if angle is 0 or 180

	inline T getAngle() const
	{
		return acos(((*this)(0, 0) + (*this)(1, 1) + (*this)(2, 2) - T(1)) / T(2));
	}

	inline T getAngleDegrees() const
	{
		return getAngle() * T(LINALG_RAD2DEG);
	}

	vec3 getAxis() const
	{
		const T m01 = (*this)[0][1];
		const T m02 = (*this)[0][2];

		const T m10 = (*this)[1][0];
		const T m20 = (*this)[2][0];

		const T m12 = (*this)[1][2];
		const T m21 = (*this)[2][1];

		const T sq = sqrtf(
			(m21 - m12) * (m21 - m12) +
			(m02 - m20) * (m02 - m20) +
			(m10 - m01) * (m10 - m01)
		);

		return vec3(
			(m21 - m12) / sq,
			(m02 - m20) / sq,
			(m10 - m01) / sq
		);
	}


	inline void swap(mat3 &other)
	{
		const mat3 tmp(*this);
		(*this) = other;
		other = tmp;
	}
	friend inline void swap(mat3 &a, mat3 &b) { a.swap(b); }
};


template<typename T>
class mat4_t
{
private:

	typedef vec2_t<T> vec2;
	typedef vec3_t<T> vec3;
	typedef vec4_t<T> vec4;

	typedef mat2_t<T> mat2;
	typedef mat3_t<T> mat3;
	typedef mat4_t<T> mat4;


public:

	static const mat4_t<T> zero;
	static const mat4_t<T> identity;


public:

	static inline mat4 translation(const vec3 &translation)
	{
		return mat4(
			vec4(T(1), T(0), T(0), T(0)),
			vec4(T(0), T(1), T(0), T(0)),
			vec4(T(0), T(0), T(1), T(0)),
			vec4(translation.x, translation.y, translation.z, T(1))
		);
	}

	static inline mat4 translation(const T tx, const T ty, const T tz = T(0))
	{
		return mat4::translation(vec3(tx, ty, tz));
	}


	static inline mat4 scaling(const vec3 &scaling)
	{
		return mat4(
			vec4(scaling.x, T(0), T(0), T(0)),
			vec4(T(0), scaling.y, T(0), T(0)),
			vec4(T(0), T(0), scaling.z, T(0)),
			vec4(T(0), T(0), T(0), T(1))
		);
	}

	static inline mat4 scaling(const T sx, const T sy, const T sz = T(1))
	{
		return mat4::scaling(vec3(sx, sy, sz));
	}


	static mat4 perspective(const T fov, const T aspect, const T zNear, const T zFar)
	{
		const T fovRad = fov * T(LINALG_DEG2RAD);

		const T range = tan(fovRad * T(0.5)) * zNear;
		const T sx = (zNear * T(2)) / (range * aspect + range * aspect);
		const T sy = zNear / range;
		const T sz = -(zFar + zNear) / (zFar - zNear);
		const T pz = -(zFar * zNear * T(2)) / (zFar - zNear);

		return mat4(
			vec4(sx, T(0), T(0), T(0)),
			vec4(T(0), sy, T(0), T(0)),
			vec4(T(0), T(0), sz, T(-1)),
			vec4(T(0), T(0), pz, T(0))
		);
	}

	static inline mat4 perspective(const T fov, const T width, const T height, const T zNear, const T zFar)
	{
		return perspective(fov, width / height, zNear, zFar);
	}

	static inline mat4 perspective(const T fov, const int width, const int height, const T zNear, const T zFar)
	{
		return perspective(fov, static_cast<T>(width), static_cast<T>(height), zNear, zFar);
	}


	static inline mat4 orthographic(const T left, const T right, const T bottom, const T top, const T zNear = T(-1), const T zFar = T(1))
	{
		return mat4(
			vec4((T(2) / (right - left)), T(0), T(0), T(0)),
			vec4(T(0), (T(2) / (top - bottom)), T(0), T(0)),
			vec4(T(0), T(0), (T(-2) / (zFar - zNear)), T(0)),
			vec4(-((right + left) / (right - left)), -((top + bottom) / (top - bottom)), -((zFar + zNear) / (zFar - zNear)), T(1))
		);
	}

	static inline mat4 ortho(const T left, const T right, const T bottom, const T top, const T zNear = T(-1), const T zFar = T(1))
	{
		return mat4::orthographic(left, right, bottom, top, zNear, zFar);
	}

	static inline mat4 ortho2d(const T left, const T right, const T bottom, const T top)
	{
		return mat4::orthographic(left, right, bottom, top, T(-1), T(1));
	}


	static inline mat4 frustum(const T left, const T right, const T bottom, const T top, const T zNear = T(-1), const T zFar = T(1))
	{
		return mat4::orthographic(left, right, bottom, top, zNear, zFar);
	}


	static mat4 viewport(const T x, const T y, const T width, const T height)
	{
		const T halfWidth = width * T(0.5);
		const T halfHeight = height * T(0.5);

		// This is only correct when glDepthRangef(0.0f, 1.0f)
		const T zNear = T(0);
		const T zFar = T(1);

		return mat4(
			halfWidth, T(0), T(0), T(0),
			T(0), halfHeight, T(0), T(0),
			T(0), T(0), (zFar - zNear) * T(0.5), T(0),
			x + halfWidth, y + halfHeight, (zNear + zFar) * T(0.5), T(1)
		);
	}


	static mat4 lookAt(const vec3 &eye, const vec3 &at, const vec3 &up = vec3::up)
	{
		const vec3 forward = normalize(at - eye);
		const vec3 right = normalize(cross(forward, up));
		const vec3 upNew = cross(right, forward);

		return mat4(
			right.x, upNew.x, -forward.x, T(0),
			right.y, upNew.y, -forward.y, T(0),
			right.z, upNew.z, -forward.z, T(0),
			-dot(right, eye), -dot(upNew, eye), -dot(-forward, eye), T(1)
		);
	}

	static inline mat4 lookAtYZ(const vec3 &eye, vec3 at, const vec3 &up = vec3::up)
	{
		return mat4::lookAt(eye, vec3(eye.x, at.y, at.z), up);
	}

	static inline mat4 lookAtXZ(const vec3 &eye, vec3 at, const vec3 &up = vec3::up)
	{
		return mat4::lookAt(eye, vec3(at.x, eye.y, at.z), up);
	}

	static inline mat4 lookAtXY(const vec3 &eye, const vec3 &at, const vec3 &up = vec3::up)
	{
		return mat4::lookAt(eye, vec3(at.x, at.y, eye.z), up);
	}


	static mat4 lookDirection(const vec3 &direction, const vec3 &up = vec3::up)
	{
		const vec3 forward = normalize(direction);
		const vec3 right = normalize(cross(forward, up));
		const vec3 upNew = cross(right, forward);

		return mat4(
			right.x, upNew.x, -forward.x, T(0),
			right.y, upNew.y, -forward.y, T(0),
			right.z, upNew.z, -forward.z, T(0),
			T(0), T(0), T(0), T(1)
		);
	}


	static vec3 project(const vec3 &object, const mat4 &modelView, const mat4 &projection, const ivec4 &viewport)
	{
		const vec4 a = modelView * vec4(object, 1.0f);
		const vec3 ab = mat3(projection) * vec3(a);

		vec4 b = vec4(ab, -a.z);

		if (LINALG_FEQUAL(b.w, 0.0f))
			return vec3(0.0f, 0.0f, 0.0f);

		b.w = 1.0f / b.w;

		b.x = b.x * b.w;
		b.y = b.y * b.w;
		b.z = b.z * b.w;

		vec3 window;

		window.x = (b.x * 0.5f + 0.5f) * viewport.z + viewport.x;
		window.y = (b.y * 0.5f + 0.5f) * viewport.w + viewport.y;

		// This is only correct when glDepthRangef(0.0f, 1.0f)
		window.z = (1.0f + b.z) * 0.5f;

		return window;
	}

	static inline vec3 project(const vec3 &object, const mat4 &model, const mat4 &view, const mat4 &projection, const ivec4 &viewport)
	{
		return mat4::project(object, (view * model), projection, viewport);
	}


	static vec3 unproject(const vec3 &window, const mat4 &modelViewProjection, const ivec4 &viewport)
	{
		const mat4 inverseMVP = mat4(modelViewProjection).inverse();

		vec4 in = vec4(window, 1.0f);

		in.x = (in.x - static_cast<float>(viewport.x)) / static_cast<float>(viewport.z);
		in.y = (in.y - static_cast<float>(viewport.y)) / static_cast<float>(viewport.w);

		in = in * 2.0f - 1.0f;
		in.w = 1.0f;

		vec4 out = inverseMVP * in;

		if (LINALG_FEQUAL(out.w, 0.0f))
			return vec3(0.0f, 0.0f, 0.0f);

		out.w = 1.0f / out.w;

		vec3 object;
		object.x = out.x * out.w;
		object.y = out.y * out.w;
		object.z = out.z * out.w;

		return object;
	}

	static inline vec3 unproject(const vec3 &window, const mat4 &modelView, const mat4 &projection, const ivec4 &viewport)
	{
		return unproject(window, (projection * modelView), viewport);
	}

	static inline vec3 unproject(const vec3 &window, const mat4 &model, const mat4 &view, const mat4 &projection, const ivec4 &viewport)
	{
		return unproject(window, (projection * view * model), viewport);
	}


	static mat4 pickMatrix(const vec2 &center, const vec2 &size, const ivec4 &viewport)
	{
		mat4 m = mat4::identity;

		if ((size.x <= T(0)) && (size.y <= T(0)))
			return m;

		m.translate(
			(viewport.z - T(2) * (center.x - viewport.x)) / size.x,
			(viewport.w - T(2) * (center.y - viewport.y)) / size.y,
			T(0));

		m.scale(
			viewport.z / size.x,
			viewport.w / size.y,
			T(1));

		return m;
	}


	static void anglesToAxes(const vec3 &angles, vec3 &left, vec3 &up, vec3 &forward)
	{
		const T sx = sinf(angles.x);
		const T cx = cosf(angles.x);

		const T sy = sinf(angles.y);
		const T cy = cosf(angles.y);

		const T sz = sinf(angles.z);
		const T cz = cosf(angles.z);

		left.x = cy * cz;
		left.y = sx * sy * cz + cx * sz;
		left.z = -cx * sy * cz + sx * sz;

		up.x = -cy * sz;
		up.y = -sx * sy * sz + cx * cz;
		up.z = cx * sy * sz + sx * cz;

		forward.x = sy;
		forward.y = -sx * cy;
		forward.z = cx * cy;
	}


public:

	vec4 columns[4];


public:

	mat4_t(const T mainDiagonalValue = T(1))
	{
		this->columns[0] = vec4(mainDiagonalValue, T(0), T(0), T(0));
		this->columns[1] = vec4(T(0), mainDiagonalValue, T(0), T(0));
		this->columns[2] = vec4(T(0), T(0), mainDiagonalValue, T(0));
		this->columns[3] = vec4(T(0), T(0), T(0), mainDiagonalValue);
	}

	mat4_t(
		const vec4 &column1, // first column
		const vec4 &column2, // second column
		const vec4 &column3, // third column
		const vec4 &column4) // fourth column
	{
		this->columns[0] = column1;
		this->columns[1] = column2;
		this->columns[2] = column3;
		this->columns[3] = column4;
	}

	mat4_t(const vec4 columns[4])
	{
		this->columns[0] = columns[0];
		this->columns[1] = columns[1];
		this->columns[2] = columns[2];
		this->columns[3] = columns[3];
	}

	mat4_t(const T values[4 * 4])
	{
		for (int i = 0; i < 4 * 4; i++)
			(reinterpret_cast<T*>(this))[i] = values[i];
	}

	mat4_t(
		const T a, const T b, const T c, const T d, // first column
		const T e, const T f, const T g, const T h, // second column
		const T i, const T j, const T k, const T l, // third column
		const T m, const T n, const T o, const T p) // fourth column
	{
		(*this) = mat4(
			vec4(a, b, c, d),
			vec4(e, f, g, h),
			vec4(i, j, k, l),
			vec4(m, n, o, p)
		);
	}

	mat4_t(const mat2 &m)
	{
		(*this) = mat4(
			vec4(m[0], T(0), T(0)),
			vec4(m[1], T(0), T(0)),
			vec4(T(0), T(0), T(1), T(0)),
			vec4(T(0), T(0), T(0), T(1))
		);
	}

	mat4_t(const mat3 &m)
	{
		(*this) = mat4(
			vec4(m[0], T(0)),
			vec4(m[1], T(0)),
			vec4(m[2], T(0)),
			vec4(T(0), T(0), T(0), T(1))
		);
	}

	~mat4_t() {}


#pragma region Operator Overloading

#pragma region Member Access Operators

	inline vec4& operator[](const int index) { return (reinterpret_cast<vec4*>(this))[index]; }
	inline vec4 operator[](const int index) const { return ((vec4*) this)[index]; }

#pragma endregion

#pragma region Function Call Operator

	T operator()(const int row, const int column) const
	{
		return (*this)[column][row];
	}

	T& operator()(const int row, const int column)
	{
		return (*this)[column][row];
	}

#pragma endregion

#pragma region Arithmetic Operators

	mat4 operator+(const mat4 &rhs) const
	{
		mat4 result;

		for (int i = 0; i < 4; i++)
			result[i] = (*this)[i] + rhs[i];

		return result;
	}

	mat4 operator-(const mat4 &rhs) const
	{
		mat4 result;

		for (int i = 0; i < 4; i++)
			result[i] = (*this)[i] - rhs[i];

		return result;
	}

	mat4 operator*(const mat4 &rhs) const
	{
		mat4 result;

		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				result.value(i, j, row(i).dot(rhs.col(j)));

		return result;
	}

	vec4 operator*(const vec4 &rhs) const
	{
		return vec4(
			(rhs.x * (*this)[0].x) + (rhs.y * (*this)[1].x) + (rhs.z * (*this)[2].x) + (rhs.w * (*this)[3].x),
			(rhs.x * (*this)[0].y) + (rhs.y * (*this)[1].y) + (rhs.z * (*this)[2].y) + (rhs.w * (*this)[3].y),
			(rhs.x * (*this)[0].z) + (rhs.y * (*this)[1].z) + (rhs.z * (*this)[2].z) + (rhs.w * (*this)[3].z),
			(rhs.x * (*this)[0].w) + (rhs.y * (*this)[1].w) + (rhs.z * (*this)[2].w) + (rhs.w * (*this)[3].w)
		);
	}

	friend vec4 operator*(const vec4 &lhs, const mat4 &rhs)
	{
		return vec4(
			lhs.dot(rhs[0]),
			lhs.dot(rhs[1]),
			lhs.dot(rhs[2]),
			lhs.dot(rhs[3])
		);
	}

	mat4 operator*(const T &rhs) const
	{
		mat4 result;

		for (int i = 0; i < 4; i++)
			result[i] = (*this)[i] * rhs;

		return result;
	}
	friend inline mat4 operator*(const T &lhs, const mat4 &rhs) { return (rhs * lhs); }

	mat4 operator/(const T &rhs) const
	{
		mat4 result;

		for (int i = 0; i < 4; i++)
			result[i] = (*this)[i] / rhs;

		return result;
	}

#pragma endregion
#pragma region Assignment Operators

	mat4& operator+=(const mat4 &rhs) { return ((*this) = (*this) + rhs); }
	mat4& operator-=(const mat4 &rhs) { return ((*this) = (*this) - rhs); }
	mat4& operator*=(const mat4 &rhs) { return ((*this) = (*this) * rhs); }
	mat4& operator*=(const T rhs) { return ((*this) = (*this) * rhs); }
	mat4& operator/=(const T rhs) { return ((*this) = (*this) / rhs); }

	mat4& operator=(const T &rhs)
	{
		(*this) = mat4(rhs);

		return (*this);
	}

	template<typename T2>
	mat4& operator=(const T2 &rhs)
	{
		(*this) = mat4(rhs);

		return (*this);
	}

	mat4& operator=(const mat4 &rhs)
	{
		for (int i = 0; i < 4; i++)
			this->columns[i] = rhs[i];

		return (*this);
	}

	template<typename T2>
	mat4& operator=(const mat4_t<T2> &rhs)
	{
		for (int i = 0; i < 4; i++)
			this->columns[i] = vec4(rhs[i]);

		return (*this);
	}

#pragma endregion

#pragma region Comparison Operators

	bool operator==(const mat4 &rhs) const;

	inline bool operator!=(const mat4 &rhs) const { return !((*this) == rhs); }

#pragma endregion

#pragma region Cast Operators

	inline explicit operator T*() const
	{
		return reinterpret_cast<T*>(this);
	}

	template<typename T2>
	inline operator mat4_t<T2>() const
	{
		return mat4_t<T2>(
			static_cast<vec4_t<T2>>(this->columns[0]),
			static_cast<vec4_t<T2>>(this->columns[1]),
			static_cast<vec4_t<T2>>(this->columns[2]),
			static_cast<vec4_t<T2>>(this->columns[3])
		);
	}

#pragma endregion

#pragma region Stream Operators

#ifdef _IOSTREAM_

	friend inline std::ostream& operator<<(std::ostream &stream, const mat4 &rhs)
	{
		return (stream << "mat4(" << rhs[0] << "," << std::endl
					   << "     " << rhs[1] << "," << std::endl
					   << "     " << rhs[2] << "," << std::endl
					   << "     " << rhs[3] << ")");
	}

	friend inline std::wostream& operator<<(std::wostream &stream, const mat4 &rhs)
	{
		return (stream << L"mat4(" << rhs[0] << L"," << std::endl
					   << L"     " << rhs[1] << L"," << std::endl
					   << L"     " << rhs[2] << L"," << std::endl
					   << L"     " << rhs[3] << L")");
	}

	friend inline std::istream& operator>>(std::istream &stream, mat4 &rhs)
	{
		rhs = mat4::identity;

		char c;

		_LINALG_IN_FIND_BEGIN();
		stream >> rhs[0];

		_LINALG_IN_FIND_NEXT();
		stream >> rhs[1];

		_LINALG_IN_FIND_NEXT();
		stream >> rhs[2];

		_LINALG_IN_FIND_NEXT();
		stream >> rhs[3];

		_LINALG_IN_FIND_END();

		return stream;
	}

	friend inline std::wistream& operator>>(std::wistream &stream, mat4 &rhs)
	{
		rhs = mat4::identity;

		wchar_t c;

		_LINALG_IN_FIND_BEGIN();
		stream >> rhs[0];

		_LINALG_IN_FIND_NEXT();
		stream >> rhs[1];

		_LINALG_IN_FIND_NEXT();
		stream >> rhs[2];

		_LINALG_IN_FIND_NEXT();
		stream >> rhs[3];

		_LINALG_IN_FIND_END();

		return stream;
	}

#endif

#pragma endregion

#pragma endregion


	T determinant() const
	{
		return (
			((*this)[0][0] * (*this)[1][1] - (*this)[0][1] * (*this)[1][0]) *
			((*this)[2][2] * (*this)[3][3] - (*this)[2][3] * (*this)[3][2]) -
			((*this)[0][0] * (*this)[1][2] - (*this)[0][2] * (*this)[1][0]) *
			((*this)[2][1] * (*this)[3][3] - (*this)[2][3] * (*this)[3][1]) +
			((*this)[0][0] * (*this)[1][3] - (*this)[0][3] * (*this)[1][0]) *
			((*this)[2][1] * (*this)[3][2] - (*this)[2][2] * (*this)[3][1]) +
			((*this)[0][1] * (*this)[1][2] - (*this)[0][2] * (*this)[1][1]) *
			((*this)[2][0] * (*this)[3][3] - (*this)[2][3] * (*this)[3][0]) -
			((*this)[0][1] * (*this)[1][3] - (*this)[0][3] * (*this)[1][1]) *
			((*this)[2][0] * (*this)[3][2] - (*this)[2][2] * (*this)[3][0]) +
			((*this)[0][2] * (*this)[1][3] - (*this)[0][3] * (*this)[1][2]) *
			((*this)[2][0] * (*this)[3][1] - (*this)[2][1] * (*this)[3][0])
		);
	}
	friend inline T determinant(const mat4 &m) { return mat4(m).determinant(); }


	mat4& inverse()
	{
		T *m = reinterpret_cast<T*>(this);
		T inv[4 * 4];

		inv[0] = m[5] * m[10] * m[15] -
			m[5] * m[11] * m[14] -
			m[9] * m[6] * m[15] +
			m[9] * m[7] * m[14] +
			m[13] * m[6] * m[11] -
			m[13] * m[7] * m[10];

		inv[4] = -m[4] * m[10] * m[15] +
			m[4] * m[11] * m[14] +
			m[8] * m[6] * m[15] -
			m[8] * m[7] * m[14] -
			m[12] * m[6] * m[11] +
			m[12] * m[7] * m[10];

		inv[8] = m[4] * m[9] * m[15] -
			m[4] * m[11] * m[13] -
			m[8] * m[5] * m[15] +
			m[8] * m[7] * m[13] +
			m[12] * m[5] * m[11] -
			m[12] * m[7] * m[9];

		inv[12] = -m[4] * m[9] * m[14] +
			m[4] * m[10] * m[13] +
			m[8] * m[5] * m[14] -
			m[8] * m[6] * m[13] -
			m[12] * m[5] * m[10] +
			m[12] * m[6] * m[9];

		inv[1] = -m[1] * m[10] * m[15] +
			m[1] * m[11] * m[14] +
			m[9] * m[2] * m[15] -
			m[9] * m[3] * m[14] -
			m[13] * m[2] * m[11] +
			m[13] * m[3] * m[10];

		inv[5] = m[0] * m[10] * m[15] -
			m[0] * m[11] * m[14] -
			m[8] * m[2] * m[15] +
			m[8] * m[3] * m[14] +
			m[12] * m[2] * m[11] -
			m[12] * m[3] * m[10];

		inv[9] = -m[0] * m[9] * m[15] +
			m[0] * m[11] * m[13] +
			m[8] * m[1] * m[15] -
			m[8] * m[3] * m[13] -
			m[12] * m[1] * m[11] +
			m[12] * m[3] * m[9];

		inv[13] = m[0] * m[9] * m[14] -
			m[0] * m[10] * m[13] -
			m[8] * m[1] * m[14] +
			m[8] * m[2] * m[13] +
			m[12] * m[1] * m[10] -
			m[12] * m[2] * m[9];

		inv[2] = m[1] * m[6] * m[15] -
			m[1] * m[7] * m[14] -
			m[5] * m[2] * m[15] +
			m[5] * m[3] * m[14] +
			m[13] * m[2] * m[7] -
			m[13] * m[3] * m[6];

		inv[6] = -m[0] * m[6] * m[15] +
			m[0] * m[7] * m[14] +
			m[4] * m[2] * m[15] -
			m[4] * m[3] * m[14] -
			m[12] * m[2] * m[7] +
			m[12] * m[3] * m[6];

		inv[10] = m[0] * m[5] * m[15] -
			m[0] * m[7] * m[13] -
			m[4] * m[1] * m[15] +
			m[4] * m[3] * m[13] +
			m[12] * m[1] * m[7] -
			m[12] * m[3] * m[5];

		inv[14] = -m[0] * m[5] * m[14] +
			m[0] * m[6] * m[13] +
			m[4] * m[1] * m[14] -
			m[4] * m[2] * m[13] -
			m[12] * m[1] * m[6] +
			m[12] * m[2] * m[5];

		inv[3] = -m[1] * m[6] * m[11] +
			m[1] * m[7] * m[10] +
			m[5] * m[2] * m[11] -
			m[5] * m[3] * m[10] -
			m[9] * m[2] * m[7] +
			m[9] * m[3] * m[6];

		inv[7] = m[0] * m[6] * m[11] -
			m[0] * m[7] * m[10] -
			m[4] * m[2] * m[11] +
			m[4] * m[3] * m[10] +
			m[8] * m[2] * m[7] -
			m[8] * m[3] * m[6];

		inv[11] = -m[0] * m[5] * m[11] +
			m[0] * m[7] * m[9] +
			m[4] * m[1] * m[11] -
			m[4] * m[3] * m[9] -
			m[8] * m[1] * m[7] +
			m[8] * m[3] * m[5];

		inv[15] = m[0] * m[5] * m[10] -
			m[0] * m[6] * m[9] -
			m[4] * m[1] * m[10] +
			m[4] * m[2] * m[9] +
			m[8] * m[1] * m[6] -
			m[8] * m[2] * m[5];

		T det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

		if (LINALG_FEQUAL(det, 0.0))
		{
			(*this) = mat4::identity;
		}
		else
		{
			det = T(1) / det;

			for (int i = 0; i < 4 * 4; i++)
				m[i] = inv[i] * det;
		}

		return (*this);
	}
	friend inline mat4 inverse(const mat4 &m) { return mat4(m).inverse(); }


	mat4& transpose()
	{
		return ((*this) = mat4(
			vec4((*this)[0].x, (*this)[1].x, (*this)[2].x, (*this)[3].x),
			vec4((*this)[0].y, (*this)[1].y, (*this)[2].y, (*this)[3].y),
			vec4((*this)[0].z, (*this)[1].z, (*this)[2].z, (*this)[3].z),
			vec4((*this)[0].w, (*this)[1].w, (*this)[2].w, (*this)[3].w)
		));
	}
	friend inline mat4 transpose(const mat4 &m) { return mat4(m).transpose(); }


	inline vec4 col(const int index) const
	{
		return (*this)[index];
	}

	inline vec4 row(const int index) const
	{
		return vec4(
			(*this)[0][index],
			(*this)[1][index],
			(*this)[2][index],
			(*this)[3][index]
		);
	}


	inline T value(const int row, const int column) const
	{
		return (*this)[column][row];
	}

	inline void value(const int row, const int column, const T &value)
	{
		(*this)[column][row] = value;
	}


	inline mat4& setZero()
	{
		return ((*this) = mat4::zero);
	}

	inline mat4& setIdentity()
	{
		return ((*this) = mat4::identity);
	}


	mat4& translate(const vec3 &translation)
	{
		return ((*this) *= mat4::translation(translation));
	}
	friend inline mat4 translate(const mat4 &m, const vec3 &translation) { return mat4(m).translate(translation); }

	inline mat4& translate(const T tx, const T ty, const T tz = T(0)) { return (*this).translate(vec3(tx, ty, tz)); }
	friend inline mat4 translate(const mat4 &m, const T tx, const T ty, const T tz = T(0)) { return mat4(m).translate(tx, ty, tz); }


	mat4& scale(const vec3 &scaling)
	{
		return ((*this) *= mat4::scaling(scaling));
	}
	friend inline mat4 scale(const mat4 &m, const vec3 &scaling) { return mat4(m).scale(scaling); }

	inline mat4& scale(const T sx, const T sy, const T sz = T(1)) { return (*this).scale(vec3(sx, sy, sz)); }
	friend inline mat4 scale(const mat4 &m, const T sx, const T sy, const T sz = T(1)) { return mat4(m).scale(sx, sy, sz); }

	inline mat4& scale(const T scaling) { return (*this).scale(vec3(scaling, scaling, scaling)); }
	friend inline mat4 scale(const mat4 &m, const T scaling) { return mat4(m).scale(scaling, scaling, scaling); }


	mat4& rotate(const T radians, const vec3 &axis)
	{
		vec3 axisNormalized = axis;

		if (!axisNormalized.isUnitVector())
			axisNormalized = normalize(axisNormalized);

		const T s = sin(radians);
		const T c = cos(radians);
		const T oc = T(1) - c;

		const T x = axisNormalized.x;
		const T y = axisNormalized.y;
		const T z = axisNormalized.z;

		return ((*this) *= mat4(
			(x * x * oc + c),
			(x * y * oc - z * s),
			(x * z * oc + y * s),
			T(0),

			(y * x * oc + z * s),
			(y * y * oc + c),
			(y * z * oc - x * s),
			T(0),

			(x * z * oc - y * s),
			(y * z * oc + x * s),
			(z * z * oc + c),
			T(0),

			T(0), T(0), T(0), T(1)
		));
	}
	friend inline mat4 rotate(const mat4 &m, const T radians, const vec3 &axis) { return mat4(m).rotate(radians, axis); }

	inline mat4& rotate(const T radians, const T ax, const T ay, const T az) { return (*this).rotate(radians, vec3(ax, ay, az)); }
	friend inline mat4 rotate(const mat4 &m, const T radians, const T ax, const T ay, const T az) { return mat4(m).rotate(radians, ax, ay, az); }


	inline mat4& rotateDegrees(const T degrees, const vec3 &axis) { return rotate(degrees * T(LINALG_DEG2RAD), axis); }
	friend inline mat4 rotateDegrees(const mat4 &m, const T degrees, const vec3 &axis) { return mat4(m).rotateDegrees(degrees, axis); }

	inline mat4& rotateDegrees(const T degrees, const T ax, const T ay, const T az) { return (*this).rotateDegrees(degrees, vec3(ax, ay, az)); }
	friend inline mat4 rotateDegrees(const mat4 &m, const T degrees, const T ax, const T ay, const T az) { return mat4(m).rotateDegrees(degrees, ax, ay, az); }


	mat4& rotateX(const T radians)
	{
		const T s = sin(radians);
		const T c = cos(radians);

		return ((*this) *= mat4(
			T(1), T(0), T(0), T(0),
			T(0), c, -s, T(0),
			T(0), s, c, T(0),
			T(0), T(0), T(0), T(1)
		));
	}
	friend inline mat4 rotateX(const mat4 &m, const T radians) { return mat4(m).rotateX(radians); }

	inline mat4& rotateXDegrees(const T degrees) { return rotateX(degrees * T(LINALG_DEG2RAD)); }
	friend inline mat4 rotateXDegrees(const mat4 &m, const T degrees) { return mat4(m).rotateXDegrees(degrees); }


	mat4& rotateY(const T radians)
	{
		const T s = sin(radians);
		const T c = cos(radians);

		return ((*this) *= mat4(
			c, T(0), s, T(0),
			T(0), T(1), T(0), T(0),
			-s, T(0), c, T(0),
			T(0), T(0), T(0), T(1)
		));
	}
	friend inline mat4 rotateY(const mat4 &m, const T radians) { return mat4(m).rotateY(radians); }

	inline mat4& rotateYDegrees(const T degrees) { return rotateY(degrees * T(LINALG_DEG2RAD)); }
	friend inline mat4 rotateYDegrees(const mat4 &m, const T degrees) { return mat4(m).rotateYDegrees(degrees); }


	mat4& rotateZ(const T radians)
	{
		const T s = sin(radians);
		const T c = cos(radians);

		return ((*this) *= mat4(
			c, -s, T(0), T(0),
			s, c, T(0), T(0),
			T(0), T(0), T(1), T(0),
			T(0), T(0), T(0), T(1)
		));
	}
	friend inline mat4 rotateZ(const mat4 &m, const T radians) { return mat4(m).rotateZ(radians); }

	inline mat4& rotateZDegrees(const T degrees) { return rotateZ(degrees * T(LINALG_DEG2RAD)); }
	friend inline mat4 rotateZDegrees(const mat4 &m, const T degrees) { return mat4(m).rotateZDegrees(degrees); }


	// Skew along the x-axis and y-axis
	mat4& skew(const T x, const T y)
	{
		return ((*this) *= mat4(
			T(1), tan(x), T(0), T(0),
			tan(y), T(1), T(0), T(0),
			T(0), T(0), T(1), T(0),
			T(0), T(0), T(0), T(1)
		));
	}
	friend inline mat4 skew(const mat4 &m, const T x, const T y) { return mat4(m).skew(x, y); }

	inline mat4& skewDegrees(const T x, const T y) { return skew(x * T(LINALG_DEG2RAD), y * T(LINALG_DEG2RAD)); }
	friend inline mat4 skewDegrees(const mat4 &m, const T x, const T y) { return mat4(m).skewDegrees(x, y); }


	// Skew along the x-axis
	mat4& skewX(const T radians)
	{
		return ((*this) *= mat4(
			T(1), tan(radians), T(0), T(0),
			T(0), T(1), T(0), T(0),
			T(0), T(0), T(1), T(0),
			T(0), T(0), T(0), T(1)
		));
	}
	friend inline mat4 skewX(const mat4 &m, const T radians) { return mat4(m).skewX(radians); }

	inline mat4& skewXDegrees(const T degrees) { return skewX(degrees * T(LINALG_DEG2RAD)); }
	friend inline mat4 skewXDegrees(const mat4 &m, const T degrees) { return mat4(m).skewXDegrees(degrees); }


	// Skew along the y-axis
	mat4& skewY(const T radians)
	{
		return ((*this) *= mat4(
			T(1), T(0), T(0), T(0),
			tan(radians), T(1), T(0), T(0),
			T(0), T(0), T(1), T(0),
			T(0), T(0), T(0), T(1)
		));
	}
	friend inline mat4 skewY(const mat4 &m, const T radians) { return mat4(m).skewY(radians); }

	inline mat4& skewYDegrees(const T degrees) { return skewY(degrees * T(LINALG_DEG2RAD)); }
	friend inline mat4 skewYDegrees(const mat4 &m, const T degrees) { return mat4(m).skewYDegrees(degrees); }


	inline mat4& setTranslation(const vec3 &translation)
	{
		return ((*this) = mat4::translation(translation));
	}

	inline mat4 setTranslation(const T tx, const T ty, const T tz = T(0))
	{
		return ((*this) = mat4::translation(vec3(tx, ty, tz)));
	}

	inline vec3 getTranslation() const
	{
		const vec4 translation = (*this)[3];

		return vec3(translation.x, translation.y, translation.z);
	}


	inline mat4& setScaling(const vec3 &scaling)
	{
		return ((*this) = mat4::scaling(scaling));
	}

	inline mat4 setScaling(const T sx, const T sy, const T sz = T(0))
	{
		return ((*this) = mat4::scaling(vec3(sx, sy, sz)));
	}

	inline vec3 getScaling() const
	{
		return vec3((*this)(0, 0), (*this)(1, 1), (*this)(2, 2));
	}


	// TODO: Gives incorrect results if angle is 0 or 180

	inline T getAngle() const
	{
		return acos(((*this)(0, 0) + (*this)(1, 1) + (*this)(2, 2) - T(1)) / T(2));
	}

	inline T getAngleDegrees() const
	{
		return getAngle() * T(LINALG_RAD2DEG);
	}

	vec3 getAxis() const
	{
		const T m01 = (*this)[0][1];
		const T m02 = (*this)[0][2];

		const T m10 = (*this)[1][0];
		const T m20 = (*this)[2][0];

		const T m12 = (*this)[1][2];
		const T m21 = (*this)[2][1];

		const T sq = sqrtf(
			(m21 - m12) * (m21 - m12) +
			(m02 - m20) * (m02 - m20) +
			(m10 - m01) * (m10 - m01)
		);

		return vec3(
			(m21 - m12) / sq,
			(m02 - m20) / sq,
			(m10 - m01) / sq
		);
	}


	inline void swap(mat4 &other)
	{
		const mat4 tmp(*this);
		(*this) = other;
		other = tmp;
	}
	friend inline void swap(mat4 &a, mat4 &b) { a.swap(b); }
};


template<typename T>
class quat_t
{
private:

	typedef vec3_t<T> vec3;
	typedef vec4_t<T> vec4;

	typedef quat_t<T> quat;


public:

	static const quat_t<T> zero;
	static const quat_t<T> one;

	static const quat_t<T> identity;


public:

	T x, y, z, w;


public:

	quat_t() : x(T(0)), y(T(0)), z(T(0)), w(T(0)) {}

	quat_t(const quat_t<T> &q) : x(q.x), y(q.y), z(q.z), w(q.w) {}
	template<typename T2> quat_t(const quat_t<T2> &q) : x(T(q.x)), y(T(q.y)), z(T(q.z)), w(T(q.w)) {}

	template<typename T2> quat_t(const vec4 &v) : x(T(v[0])), y(T(v[1])), z(T(v[2])), w(T(v[3])) {}
	template<typename T2> quat_t(const vec4_t<T2> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(T(v.w)) {}

	template<typename T2> quat_t(const T2 &xyzw) : x(T(xyzw)), y(T(xyzw)), z(T(xyzw)), w(T(xyzw)) {}
	template<typename T2> quat_t(const T2 &x, const T2 &y, const T2 &z, const T2 &w) : x(T(x)), y(T(y)), z(T(z)), w(T(w)) {}

	template<typename T2> quat_t(const T2 *xyzw) : x(T(xyzw[0])), y(T(xyzw[1])), z(T(xyzw[2])), w(T(xyzw[3])) {}

	quat_t(const T radians, const vec3 &axis)
	{
		(*this) = rotate(radians, axis);
	}

	~quat_t() {}


#pragma region Operator Overloading

#pragma region Member Access Operators

	inline T& operator[](const int index) { return (reinterpret_cast<T*>(this))[index]; }
	inline T operator[](const int index) const { return ((T*) this)[index]; }

#pragma endregion

#pragma region Arithmetic Operators

	quat operator+(const quat &rhs) const { return quat(this->x + rhs.x, this->y + rhs.y, this->z + rhs.z, this->w + rhs.w); }
	quat operator-(const quat &rhs) const { return quat(this->x - rhs.x, this->y - rhs.y, this->z - rhs.z, this->w - rhs.w); }

	// Like for matrix multiplication, quaternion multiplication is non-commutative:
	// (q1 * q2) != (q2 * q1)
	quat operator*(const quat &rhs) const
	{
		return quat(
			(this->w * rhs.x) + (this->x * rhs.w) + (this->y * rhs.z) - (this->z * rhs.y),
			(this->w * rhs.y) - (this->x * rhs.z) + (this->y * rhs.w) + (this->z * rhs.x),
			(this->w * rhs.z) + (this->x * rhs.y) - (this->y * rhs.x) + (this->z * rhs.w),
			(this->w * rhs.w) - (this->x * rhs.x) - (this->y * rhs.y) - (this->z * rhs.z)
		);
	}

	friend inline quat& operator*(const quat &lhs, const T &rhs) { return (lhs * quat(rhs)); }
	friend inline quat& operator/(const quat &lhs, const T &rhs) { return (lhs / quat(rhs)); }

#pragma endregion
#pragma region Assignment Operators

	inline quat& operator+=(const quat &rhs) { return ((*this) = ((*this) + rhs)); }
	inline quat& operator-=(const quat &rhs) { return ((*this) = ((*this) - rhs)); }
	inline quat& operator*=(const quat &rhs) { return ((*this) = ((*this) * rhs)); }
	inline quat& operator*=(const T &rhs) { return ((*this) = ((*this) * rhs)); }
	inline quat& operator/=(const T &rhs) { return ((*this) = ((*this) / rhs)); }

	quat& operator=(const quat &rhs)
	{
		this->x = rhs.x;
		this->y = rhs.y;
		this->z = rhs.z;
		this->w = rhs.w;

		return (*this);
	}

	template<typename T2>
	vec4& operator=(const quat_t<T2> &rhs)
	{
		this->x = T(rhs.x);
		this->y = T(rhs.y);
		this->z = T(rhs.z);
		this->w = T(rhs.w);

		return (*this);
	}

#pragma endregion

#pragma region Cast Operators

	inline explicit operator T*() const
	{
		return reinterpret_cast<T*>(this);
	}

	template<typename T2>
	inline operator quat_t<T2>() const
	{
		return quat_t<T2>(
			static_cast<T2>(this->x),
			static_cast<T2>(this->y),
			static_cast<T2>(this->z),
			static_cast<T2>(this->w)
		);
	}

	template<typename T2>
	inline operator vec4_t<T2>() const
	{
		return vec4_t<T2>(
			static_cast<T2>(this->x),
			static_cast<T2>(this->y),
			static_cast<T2>(this->z),
			static_cast<T2>(this->w)
		);
	}

#pragma endregion

#pragma region Stream Operators

#ifdef _IOSTREAM_

	friend inline std::ostream& operator<<(std::ostream &stream, const quat &rhs)
	{
		return (stream << "quat(" << rhs.x << ", " << rhs.y << ", " << rhs.z << ", " << rhs.w << ")");
	}

	friend inline std::wostream& operator<<(std::wostream &stream, const quat &rhs)
	{
		return (stream << L"quat(" << rhs.x << L", " << rhs.y << L", " << rhs.z << L", " << rhs.w << L")");
	}

	friend inline std::istream& operator>>(std::istream &stream, quat &rhs)
	{
		rhs = quat::zero;

		char c;

		_LINALG_IN_FIND_BEGIN();
		stream >> rhs.x;

		_LINALG_IN_FIND_NEXT();
		stream >> rhs.y;

		_LINALG_IN_FIND_NEXT();
		stream >> rhs.z;

		_LINALG_IN_FIND_NEXT();
		stream >> rhs.w;

		_LINALG_IN_FIND_END();

		return stream;
	}

	friend inline std::wistream& operator>>(std::wistream &stream, quat &rhs)
	{
		rhs = quat::zero;

		wchar_t c;

		_LINALG_IN_FIND_BEGIN();
		stream >> rhs.x;

		_LINALG_IN_FIND_NEXT();
		stream >> rhs.y;

		_LINALG_IN_FIND_NEXT();
		stream >> rhs.z;

		_LINALG_IN_FIND_NEXT();
		stream >> rhs.w;

		_LINALG_IN_FIND_END();

		return stream;
	}

#endif

#pragma endregion

#pragma endregion


	quat conjugate() const
	{
		return quat(
			-this->x,
			-this->y,
			-this->z,
			this->w
		);
	}
	friend inline quat conjugate(const quat &q) { return q.conjugate(); }


	quat normalize() const
	{
		T norm = sqrt(this->x * this->x + this->y * this->y + this->z * this->z + this->w * this->w);
		norm = T(1.0) / norm;

		return quat(
			this->x * norm,
			this->y * norm,
			this->z * norm,
			this->w * norm
		);
	}
	friend inline quat normalize(const quat &q) { return q.normalize(); }


	quat rotate(const T radians, const vec3 axis) const
	{
		const T half_angle = radians * T(0.5);

		const T s = sin(half_angle);
		const T c = cos(half_angle);

		return quat(
			axis.x * s,
			axis.y * s,
			axis.z * s,
			c
		);
	}
	friend inline quat rotate(const quat &q, const T angle, const vec3 axis) { return q.rotate(angle, axis); }

	inline quat& rotateDegrees(const T degrees) { return rotate(degrees * T(LINALG_DEG2RAD)); }
	friend inline quat rotateDegrees(const quat &q, const T degrees) { return quat(q).rotateDegrees(degrees); }


	// Normalize then conjugate the current quaternion,
	// to get the inverse quaternion.
	inline quat inverse() const
	{
		return normalize().conjugate();
	}
	friend inline quat inverse(const quat &q) { return q.inverse(); }


	inline quat slerp(const quat &to, const quat &t) const
	{
		const T EPSILON = T(1E-6f);

		const quat a = (*this);
		const quat b = to;

		T omega = T(0);
		T cosom = (a.x * b.x) + (a.y * b.y) + (a.z * b.z) + (a.w * b.w);
		T sinom = T(0);
		T scale0 = T(0);
		T scale1 = T(0);

		quat result;

		if ((T(1) + cosom) > EPSILON)
		{
			// a and b quaternions are not opposite each other
			if ((T(1) - cosom) > EPSILON)
			{
				// Standard case - slerp
				omega = acosf(cosom);
				sinom = sinf(omega);
				scale0 = sinf((T(1) - t) * omega) / sinom;
				scale1 = sinf(t * omega) / sinom;
			}
			else
			{
				// a and b quaternions are very close so lerp instead
				scale0 = T(1) - t;
				scale1 = t;
			}

			result.x = scale0 * a.x + scale1 * b.x;
			result.y = scale0 * a.y + scale1 * b.y;
			result.z = scale0 * a.z + scale1 * b.z;
			result.w = scale0 * a.w + scale1 * b.w;
		}
		else
		{
			// a and b quaternions are opposite each other
			result.x = -b.y;
			result.y = b.x;
			result.z = -b.w;
			result.w = b.z;

			scale0 = sinf((T(1) - t) - T(1.57079632679));
			scale1 = sinf(t * T(1.57079632679));

			result.x = scale0 * a.x + scale1 * result.x;
			result.y = scale0 * a.y + scale1 * result.y;
			result.z = scale0 * a.z + scale1 * result.z;
			result.w = scale0 * a.w + scale1 * result.w;
		}

		return result;
	}
	friend inline quat slerp(const quat &from, const quat &to, const T t) { return from.lerp(to, t); }


	inline T getAngle() const
	{
		return T(2) * acos(this->w);
	}

	inline T getAngleDegrees() const
	{
		return getAngle() * T(LINALG_RAD2DEG);
	}

	vec3 getAxis() const
	{
		const T squared = T(1) - this->w * this->w;

		if (LINALG_FEQUAL(squared, 0.0f))
			return vec3(T(1), T(0), T(0));

		const T length = sqrt(squared);
		const T invLength = T(1) / length;

		return vec3(this->x * invLength, this->y * invLength, this->z * invLength);
	}


	inline void swap(quat &other)
	{
		const quat tmp(*this);
		(*this) = other;
		other = tmp;
	}
	friend inline void swap(quat &a, quat &b) { a.swap(b); }
};


// It isn't an optimal solution, to inline all template functions that has explicit specialization.
// But it is needed if we don't want to run into "multiple definitions" compilation error.


#pragma region vec2

#pragma region Static Members

template<typename T> const vec2_t<T> vec2_t<T>::zero = vec2_t<T>(T(0), T(0));
template<typename T> const vec2_t<T> vec2_t<T>::one = vec2_t<T>(T(1), T(1));

template<typename T> const vec2_t<T> vec2_t<T>::up = vec2_t<T>(T(0), T(1));
template<typename T> const vec2_t<T> vec2_t<T>::down = vec2_t<T>(T(0), T(-1));

template<typename T> const vec2_t<T> vec2_t<T>::left = vec2_t<T>(T(-1), T(0));
template<typename T> const vec2_t<T> vec2_t<T>::right = vec2_t<T>(T(1), T(0));

#pragma endregion

#pragma region Constructors

template<typename T>
template<typename T2>
vec2_t<T>::vec2_t(const vec3_t<T2> &v) : x(T(v.x)), y(T(v.y))
{

}

template<typename T>
template<typename T2>
vec2_t<T>::vec2_t(const vec4_t<T2> &v) : x(T(v.x)), y(T(v.y))
{

}

#pragma endregion

#pragma region Comparison Operators

template<typename T> inline bool vec2_t<T>::operator==(const vec2 &rhs) const
{
	return ((this->x == rhs.x) && (this->y == rhs.y));
}

template<> inline bool fvec2::operator==(const fvec2 &rhs) const
{
	return (LINALG_FEQUAL(this->x, rhs.x) && LINALG_FEQUAL(this->y, rhs.y));
}

template<> inline bool dvec2::operator==(const dvec2 &rhs) const
{
	return (LINALG_DEQUAL(this->x, rhs.x) && LINALG_DEQUAL(this->y, rhs.y));
}

#pragma endregion

#pragma region Normalize

template<> inline fvec2 fvec2::normalize(const float &to) const
{
	float len = length();

	if (!LINALG_FEQUAL(len, 0.0f) && !LINALG_FEQUAL(len, to))
	{
		len = to / len;

		return fvec2(*this) * len;
	}

	return (*this);
}

template<> inline dvec2 dvec2::normalize(const double &to) const
{
	double len = length();

	if (!LINALG_DEQUAL(len, 0.0) && !LINALG_DEQUAL(len, to))
	{
		len = to / len;

		return dvec2(*this) * len;
	}

	return (*this);
}

#pragma endregion

#pragma region Check Vector Kind

template<typename T> bool vec2_t<T>::isNullVector() const
{
	// It's faster to just compare to 0, than to calculate the length
	return ((this->x == 0) && (this->y == 0));
}

template<> bool fvec2::isNullVector() const
{
	// It's faster to just compare to 0, than to calculate the length
	return (LINALG_FEQUAL(this->x, 0.0f) && LINALG_FEQUAL(this->y, 0.0f));
}

template<> bool dvec2::isNullVector() const
{
	// It's faster to just compare to 0, than to calculate the length
	return (LINALG_DEQUAL(this->x, 0.0) && LINALG_DEQUAL(this->y, 0.0));
}


template<> inline bool fvec2::isUnitVector() const
{
	return LINALG_FEQUAL(lengthSquared(), 1.0f);
}

template<> inline bool dvec2::isUnitVector() const
{
	return LINALG_DEQUAL(lengthSquared(), 1.0f);
}


template<> bool fvec2::isNormalized(const float &to) const
{
	const float len = length();

	return LINALG_FEQUAL(len, to);
}

template<> bool dvec2::isNormalized(const double &to) const
{
	const double len = length();

	return LINALG_DEQUAL(len, to);
}


template<typename T> bool vec2_t<T>::isOrthogonalTo(const vec2 &rhs) const
{
	return (dot(rhs) == 0);
}

template<> bool fvec2::isOrthogonalTo(const vec2 &rhs) const
{
	const float d = dot(rhs);

	return (LINALG_FEQUAL(d, 0.0f));
}

template<> bool dvec2::isOrthogonalTo(const vec2 &rhs) const
{
	const double d = dot(rhs);

	return (LINALG_DEQUAL(d, 0.0));
}


template<typename T> bool vec2_t<T>::isPerpendicularTo(const vec2 &rhs) const
{
	return (dot(rhs) == 0);
}

template<> bool fvec2::isPerpendicularTo(const vec2 &rhs) const
{
	const float d = dot(rhs);

	return (LINALG_FEQUAL(d, 0.0f));
}

template<> bool dvec2::isPerpendicularTo(const vec2 &rhs) const
{
	const double d = dot(rhs);

	return (LINALG_DEQUAL(d, 0.0));
}


template<typename T> bool vec2_t<T>::isParallelTo(const vec2 &rhs) const
{
	return (dot(rhs) == 1);
}

template<> bool fvec2::isParallelTo(const vec2 &rhs) const
{
	const float d = dot(rhs);

	return (LINALG_FEQUAL(d, 1.0f));
}

template<> bool dvec2::isParallelTo(const vec2 &rhs) const
{
	const double d = dot(rhs);

	return (LINALG_DEQUAL(d, 1.0));
}

#pragma endregion

#pragma region Validate sizeof Templated Objects

#ifdef DEBUG

#ifndef STATIC_ASSERT
#	define STATIC_ASSERT(bool_constexpr) static_assert(bool_constexpr, #bool_constexpr)
#endif

STATIC_ASSERT(sizeof(vec2_t<float>) == (sizeof(float) * 2));
STATIC_ASSERT(sizeof(vec2_t<double>) == (sizeof(double) * 2));

STATIC_ASSERT(sizeof(vec2_t<signed int>) == (sizeof(signed int) * 2));
STATIC_ASSERT(sizeof(vec2_t<unsigned int>) == (sizeof(unsigned int) * 2));

STATIC_ASSERT(sizeof(vec2_t<bool>) == (sizeof(bool) * 2));

STATIC_ASSERT(sizeof(vec2_t<signed long>) == (sizeof(signed long) * 2));
STATIC_ASSERT(sizeof(vec2_t<unsigned long>) == (sizeof(unsigned long) * 2));

STATIC_ASSERT(sizeof(vec2_t<signed long long>) == (sizeof(signed long long) * 2));
STATIC_ASSERT(sizeof(vec2_t<unsigned long long>) == (sizeof(unsigned long long) * 2));

#endif

#pragma endregion

#pragma endregion


#pragma region vec3

#pragma region Static Members

template<typename T> const vec3_t<T> vec3_t<T>::zero = vec3_t<T>(T(0), T(0), T(0));
template<typename T> const vec3_t<T> vec3_t<T>::one = vec3_t<T>(T(1), T(1), T(1));

template<typename T> const vec3_t<T> vec3_t<T>::up = vec3_t<T>(T(0), T(1), T(0));
template<typename T> const vec3_t<T> vec3_t<T>::down = vec3_t<T>(T(0), T(-1), T(0));

template<typename T> const vec3_t<T> vec3_t<T>::left = vec3_t<T>(T(-1), T(0), T(0));
template<typename T> const vec3_t<T> vec3_t<T>::right = vec3_t<T>(T(1), T(0), T(0));

template<typename T> const vec3_t<T> vec3_t<T>::forward = vec3_t<T>(T(0), T(0), T(1));
template<typename T> const vec3_t<T> vec3_t<T>::backward = vec3_t<T>(T(0), T(0), T(-1));

#pragma endregion

#pragma region Constructors

template<typename T>
template<typename T2>
vec3_t<T>::vec3_t(const vec4_t<T2> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z))
{

}

#pragma endregion

#pragma region Comparison Operators

template<typename T> inline bool vec3_t<T>::operator==(const vec3 &rhs) const
{
	return ((this->x == rhs.x) && (this->y == rhs.y) && (this->z == rhs.z));
}

template<> inline bool fvec3::operator==(const fvec3 &rhs) const
{
	return (LINALG_FEQUAL(this->x, rhs.x) && LINALG_FEQUAL(this->y, rhs.y) && LINALG_FEQUAL(this->z, rhs.z));
}

template<> inline bool dvec3::operator==(const dvec3 &rhs) const
{
	return (LINALG_DEQUAL(this->x, rhs.x) && LINALG_DEQUAL(this->y, rhs.y) && LINALG_DEQUAL(this->z, rhs.z));
}

#pragma endregion

#pragma region Normalize

template<> inline fvec3 fvec3::normalize(const float &to) const
{
	float len = length();

	if (!LINALG_FEQUAL(len, 0.0f) && !LINALG_FEQUAL(len, to))
	{
		len = to / len;

		return fvec3(*this) * len;
	}

	return (*this);
}

template<> inline dvec3 dvec3::normalize(const double &to) const
{
	double len = length();

	if (!LINALG_DEQUAL(len, 0.0) && !LINALG_DEQUAL(len, to))
	{
		len = to / len;

		return dvec3(*this) * len;
	}

	return (*this);
}

#pragma endregion

#pragma region Check Vector Kind

template<typename T> bool vec3_t<T>::isNullVector() const
{
	// It's faster to just compare to 0, than to calculate the length
	return ((this->x == 0) && (this->y == 0) && (this->z == 0));
}

template<> bool fvec3::isNullVector() const
{
	// It's faster to just compare to 0, than to calculate the length
	return (LINALG_FEQUAL(this->x, 0.0f) && LINALG_FEQUAL(this->y, 0.0f) && LINALG_FEQUAL(this->z, 0.0f));
}

template<> bool dvec3::isNullVector() const
{
	// It's faster to just compare to 0, than to calculate the length
	return (LINALG_DEQUAL(this->x, 0.0) && LINALG_DEQUAL(this->y, 0.0) && LINALG_DEQUAL(this->z, 0.0));
}


template<> inline bool fvec3::isUnitVector() const
{
	return LINALG_FEQUAL(lengthSquared(), 1.0f);
}

template<> inline bool dvec3::isUnitVector() const
{
	return LINALG_DEQUAL(lengthSquared(), 1.0f);
}


template<> bool fvec3::isNormalized(const float &to) const
{
	const float len = length();

	return LINALG_FEQUAL(len, to);
}

template<> bool dvec3::isNormalized(const double &to) const
{
	const double len = length();

	return LINALG_DEQUAL(len, to);
}


template<typename T> bool vec3_t<T>::isOrthogonalTo(const vec3 &rhs) const
{
	return (dot(rhs) == 0);
}

template<> bool fvec3::isOrthogonalTo(const vec3 &rhs) const
{
	const float d = dot(rhs);

	return (LINALG_FEQUAL(d, 0.0f));
}

template<> bool dvec3::isOrthogonalTo(const vec3 &rhs) const
{
	const double d = dot(rhs);

	return (LINALG_DEQUAL(d, 0.0));
}


template<typename T> bool vec3_t<T>::isPerpendicularTo(const vec3 &rhs) const
{
	return (dot(rhs) == 0);
}

template<> bool fvec3::isPerpendicularTo(const vec3 &rhs) const
{
	const float d = dot(rhs);

	return (LINALG_FEQUAL(d, 0.0f));
}

template<> bool dvec3::isPerpendicularTo(const vec3 &rhs) const
{
	const double d = dot(rhs);

	return (LINALG_DEQUAL(d, 0.0));
}


template<typename T> bool vec3_t<T>::isParallelTo(const vec3 &rhs) const
{
	return (dot(rhs) == 1);
}

template<> bool fvec3::isParallelTo(const vec3 &rhs) const
{
	const float d = dot(rhs);

	return (LINALG_FEQUAL(d, 1.0f));
}

template<> bool dvec3::isParallelTo(const vec3 &rhs) const
{
	const double d = dot(rhs);

	return (LINALG_DEQUAL(d, 1.0));
}

#pragma endregion

#pragma region Validate sizeof Templated Objects

#ifdef DEBUG

#ifndef STATIC_ASSERT
#	define STATIC_ASSERT(bool_constexpr) static_assert(bool_constexpr, #bool_constexpr)
#endif

STATIC_ASSERT(sizeof(vec3_t<float>) == (sizeof(float) * 3));
STATIC_ASSERT(sizeof(vec3_t<double>) == (sizeof(double) * 3));

STATIC_ASSERT(sizeof(vec3_t<signed int>) == (sizeof(signed int) * 3));
STATIC_ASSERT(sizeof(vec3_t<unsigned int>) == (sizeof(unsigned int) * 3));

STATIC_ASSERT(sizeof(vec3_t<bool>) == (sizeof(bool) * 3));

STATIC_ASSERT(sizeof(vec3_t<signed long>) == (sizeof(signed long) * 3));
STATIC_ASSERT(sizeof(vec3_t<unsigned long>) == (sizeof(unsigned long) * 3));

STATIC_ASSERT(sizeof(vec3_t<signed long long>) == (sizeof(signed long long) * 3));
STATIC_ASSERT(sizeof(vec3_t<unsigned long long>) == (sizeof(unsigned long long) * 3));

#endif

#pragma endregion

#pragma endregion


#pragma region vec4

#pragma region Static Members

template<typename T> const vec4_t<T> vec4_t<T>::zero = vec4_t<T>(T(0), T(0), T(0), T(0));
template<typename T> const vec4_t<T> vec4_t<T>::one = vec4_t<T>(T(1), T(1), T(1), T(1));

#pragma endregion

#pragma region Constructors



#pragma endregion

#pragma region Comparison Operators

template<typename T> inline bool vec4_t<T>::operator==(const vec4 &rhs) const
{
	return ((this->x == rhs.x) && (this->y == rhs.y) && (this->z == rhs.z) && (this->w == rhs.w));
}

template<> inline bool fvec4::operator==(const fvec4 &rhs) const
{
	return (LINALG_FEQUAL(this->x, rhs.x) && LINALG_FEQUAL(this->y, rhs.y) && LINALG_FEQUAL(this->z, rhs.z) && LINALG_FEQUAL(this->w, rhs.w));
}

template<> inline bool dvec4::operator==(const dvec4 &rhs) const
{
	return (LINALG_DEQUAL(this->x, rhs.x) && LINALG_DEQUAL(this->y, rhs.y) && LINALG_DEQUAL(this->z, rhs.z) && LINALG_DEQUAL(this->w, rhs.w));
}

#pragma endregion

#pragma region Normalize

template<> inline fvec4 fvec4::normalize(const float &to) const
{
	float len = length();

	if (!LINALG_FEQUAL(len, 0.0f) && !LINALG_FEQUAL(len, to))
	{
		len = to / len;

		return fvec4(*this) * len;
	}

	return (*this);
}

template<> inline dvec4 dvec4::normalize(const double &to) const
{
	double len = length();

	if (!LINALG_DEQUAL(len, 0.0) && !LINALG_DEQUAL(len, to))
	{
		len = to / len;

		return dvec4(*this) * len;
	}

	return (*this);
}

#pragma endregion

#pragma region Check Vector Kind

template<typename T> bool vec4_t<T>::isNullVector() const
{
	// It's faster to just compare to 0, than to calculate the length
	return ((this->x == 0) && (this->y == 0) && (this->z == 0) && (this->w == 0));
}

template<> bool fvec4::isNullVector() const
{
	// It's faster to just compare to 0, than to calculate the length
	return (LINALG_FEQUAL(this->x, 0.0f) && LINALG_FEQUAL(this->y, 0.0f) && LINALG_FEQUAL(this->z, 0.0f) && LINALG_FEQUAL(this->w, 0.0f));
}

template<> bool dvec4::isNullVector() const
{
	// It's faster to just compare to 0, than to calculate the length
	return (LINALG_DEQUAL(this->x, 0.0) && LINALG_DEQUAL(this->y, 0.0) && LINALG_DEQUAL(this->z, 0.0) && LINALG_DEQUAL(this->w, 0.0));
}


template<> inline bool fvec4::isUnitVector() const
{
	return LINALG_FEQUAL(lengthSquared(), 1.0f);
}

template<> inline bool dvec4::isUnitVector() const
{
	return LINALG_DEQUAL(lengthSquared(), 1.0f);
}


template<> bool fvec4::isNormalized(const float &to) const
{
	const float len = length();

	return LINALG_FEQUAL(len, to);
}

template<> bool dvec4::isNormalized(const double &to) const
{
	const double len = length();

	return LINALG_DEQUAL(len, to);
}


template<typename T> bool vec4_t<T>::isOrthogonalTo(const vec4 &rhs) const
{
	return (dot(rhs) == 0);
}

template<> bool fvec4::isOrthogonalTo(const vec4 &rhs) const
{
	const float d = dot(rhs);

	return (LINALG_FEQUAL(d, 0.0f));
}

template<> bool dvec4::isOrthogonalTo(const vec4 &rhs) const
{
	const double d = dot(rhs);

	return (LINALG_DEQUAL(d, 0.0));
}


template<typename T> bool vec4_t<T>::isPerpendicularTo(const vec4 &rhs) const
{
	return (dot(rhs) == 0);
}

template<> bool fvec4::isPerpendicularTo(const vec4 &rhs) const
{
	const float d = dot(rhs);

	return (LINALG_FEQUAL(d, 0.0f));
}

template<> bool dvec4::isPerpendicularTo(const vec4 &rhs) const
{
	const double d = dot(rhs);

	return (LINALG_DEQUAL(d, 0.0));
}


template<typename T> bool vec4_t<T>::isParallelTo(const vec4 &rhs) const
{
	return (dot(rhs) == 1);
}

template<> bool fvec4::isParallelTo(const vec4 &rhs) const
{
	const float d = dot(rhs);

	return (LINALG_FEQUAL(d, 1.0f));
}

template<> bool dvec4::isParallelTo(const vec4 &rhs) const
{
	const double d = dot(rhs);

	return (LINALG_DEQUAL(d, 1.0));
}

#pragma endregion

#pragma region Validate sizeof Templated Objects

#ifdef DEBUG

#ifndef STATIC_ASSERT
#	define STATIC_ASSERT(bool_constexpr) static_assert(bool_constexpr, #bool_constexpr)
#endif

STATIC_ASSERT(sizeof(vec4_t<float>) == (sizeof(float) * 4));
STATIC_ASSERT(sizeof(vec4_t<double>) == (sizeof(double) * 4));

STATIC_ASSERT(sizeof(vec4_t<signed int>) == (sizeof(signed int) * 4));
STATIC_ASSERT(sizeof(vec4_t<unsigned int>) == (sizeof(unsigned int) * 4));

STATIC_ASSERT(sizeof(vec4_t<bool>) == (sizeof(bool) * 4));

STATIC_ASSERT(sizeof(vec4_t<signed long>) == (sizeof(signed long) * 4));
STATIC_ASSERT(sizeof(vec4_t<unsigned long>) == (sizeof(unsigned long) * 4));

STATIC_ASSERT(sizeof(vec4_t<signed long long>) == (sizeof(signed long long) * 4));
STATIC_ASSERT(sizeof(vec4_t<unsigned long long>) == (sizeof(unsigned long long) * 4));

#endif

#pragma endregion

#pragma endregion


#pragma region mat2

#pragma region Static Members

template<typename T> const mat2_t<T> mat2_t<T>::zero = mat2_t<T>(T(0));
template<typename T> const mat2_t<T> mat2_t<T>::identity = mat2_t<T>(T(1));

#pragma endregion

#pragma region Constructors

template<typename T>
template<typename T2>
mat2_t<T>::mat2_t(const mat3_t<T2> &m)
{
	this->columns[0] = vec2(m.columns[0]);
	this->columns[1] = vec2(m.columns[1]);
}

template<typename T>
template<typename T2>
mat2_t<T>::mat2_t(const mat4_t<T2> &m)
{
	this->columns[0] = vec2(m.columns[0]);
	this->columns[1] = vec2(m.columns[1]);
}

#pragma endregion

#pragma region Comparison Operators

template<typename T> inline bool mat2_t<T>::operator==(const mat2_t &rhs) const
{
	for (int i = 0; i < 2; i++)
		if ((*this)[i] != rhs[i])
			return false;

	return true;
}

template<> inline bool mat2_t<float>::operator==(const mat2_t<float> &rhs) const
{
	for (int i = 0; i < 2; i++)
		if (!LINALG_FEQUAL((*this)[i], rhs[i]))
			return false;

	return true;
}

template<> inline bool mat2_t<double>::operator==(const mat2_t<double> &rhs) const
{
	for (int i = 0; i < 2; i++)
		if (!LINALG_DEQUAL((*this)[i], rhs[i]))
			return false;

	return true;
}

#pragma endregion

#pragma endregion


#pragma region mat3

#pragma region Static Members

template<typename T> const mat3_t<T> mat3_t<T>::zero = mat3_t<T>(T(0));
template<typename T> const mat3_t<T> mat3_t<T>::identity = mat3_t<T>(T(1));

#pragma endregion

#pragma region Constructors

template<typename T>
template<typename T2>
mat3_t<T>::mat3_t(const mat4_t<T2> &m)
{
	this->columns[0] = vec3(m.columns[0]);
	this->columns[1] = vec3(m.columns[1]);
	this->columns[2] = vec3(m.columns[2]);
}

#pragma endregion

#pragma region Comparison Operators

template<typename T> inline bool mat3_t<T>::operator==(const mat3_t &rhs) const
{
	for (int i = 0; i < 3; i++)
		if ((*this)[i] != rhs[i])
			return false;

	return true;
}

template<> inline bool mat3_t<float>::operator==(const mat3_t<float> &rhs) const
{
	for (int i = 0; i < 3; i++)
		if (!LINALG_FEQUAL((*this)[i], rhs[i]))
			return false;

	return true;
}

template<> inline bool mat3_t<double>::operator==(const mat3_t<double> &rhs) const
{
	for (int i = 0; i < 3; i++)
		if (!LINALG_DEQUAL((*this)[i], rhs[i]))
			return false;

	return true;
}

#pragma endregion

#pragma endregion


#pragma region mat4

#pragma region Static Members

template<typename T> const mat4_t<T> mat4_t<T>::zero = mat4_t<T>(T(0));
template<typename T> const mat4_t<T> mat4_t<T>::identity = mat4_t<T>(T(1));

#pragma endregion

#pragma region Constructors



#pragma endregion

#pragma region Comparison Operators

template<typename T> inline bool mat4_t<T>::operator==(const mat4_t &rhs) const
{
	for (int i = 0; i < 4; i++)
		if ((*this)[i] != rhs[i])
			return false;

	return true;
}

template<> inline bool mat4_t<float>::operator==(const mat4_t<float> &rhs) const
{
	for (int i = 0; i < 4; i++)
		if (!LINALG_FEQUAL((*this)[i], rhs[i]))
			return false;

	return true;
}

template<> inline bool mat4_t<double>::operator==(const mat4_t<double> &rhs) const
{
	for (int i = 0; i < 4; i++)
		if (!LINALG_DEQUAL((*this)[i], rhs[i]))
			return false;

	return true;
}

#pragma endregion

#pragma endregion


#pragma region quat

#pragma region Static Members

template<typename T> const quat_t<T> quat_t<T>::zero = vec4_t<T>(T(0), T(0), T(0), T(0));
template<typename T> const quat_t<T> quat_t<T>::one = vec4_t<T>(T(1), T(1), T(1), T(1));

template<typename T> const quat_t<T> quat_t<T>::identity = vec4_t<T>(T(0), T(0), T(0), T(1));

#pragma endregion

#pragma endregion


// Enable structure padding
#pragma pack(pop)


#endif