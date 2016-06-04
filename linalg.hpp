
// Author: Christian Vallentin <mail@vallentinsource.com>
// Website: http://vallentinsource.com
// Repository: https://github.com/MrVallentin/LinearAlgebra
//
// Date Created: October 01, 2013
// Last Modified: June 04, 2016

#ifndef LINEAR_ALGEBRA_HPP
#define LINEAR_ALGEBRA_HPP


#include <math.h>


#ifndef LINALG_DEFAULT_SCALAR
#	define LINALG_DEFAULT_SCALAR float
#endif


// Disable structure padding
#pragma pack(push, 1)


template<typename T> class vec2_t;
template<typename T> class vec3_t;
template<typename T> class vec4_t;


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


#if defined(_DEBUG) && !defined(DEBUG)
#	define DEBUG 1
#endif


// #define LINALG_FEQUAL(x, y) ((((y) - 1E-6f) < (x)) && ((x) < ((y) + 1E-6f)))
// #define LINALG_DEQUAL(x, y) ((((y) - 1E-6) < (x)) && ((x) < ((y) + 1E-6)))

// This was changed from 1E-6 to 1E-4 as asserting rotate(90deg) didn't match
#define LINALG_FEQUAL(x, y) ((((y) - 1E-4f) < (x)) && ((x) < ((y) + 1E-4f)))
#define LINALG_DEQUAL(x, y) ((((y) - 1E-4) < (x)) && ((x) < ((y) + 1E-4)))


template<typename T>
class vec2_t
{
private:

	typedef vec2_t<T> vec2;


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

	vec2_t(void) : x(T(0)), y(T(0)) {}

	// vec2_t(const vec2 &v) : x(v.x), y(v.y) {}

	vec2_t(const vec2_t<float> &v) : x(T(v.x)), y(T(v.y)) {}
	vec2_t(const vec2_t<double> &v) : x(T(v.x)), y(T(v.y)) {}
	vec2_t(const vec2_t<signed int> &v) : x(T(v.x)), y(T(v.y)) {}
	vec2_t(const vec2_t<unsigned int> &v) : x(T(v.x)), y(T(v.y)) {}
	vec2_t(const vec2_t<bool> &v) : x(T(v.x)), y(T(v.y)) {}

	vec2_t(const vec2_t<signed long> &v) : x(T(v.x)), y(T(v.y)) {}
	vec2_t(const vec2_t<unsigned long> &v) : x(T(v.x)), y(T(v.y)) {}
	vec2_t(const vec2_t<signed long long> &v) : x(T(v.x)), y(T(v.y)) {}
	vec2_t(const vec2_t<unsigned long long> &v) : x(T(v.x)), y(T(v.y)) {}

	vec2_t(const T &xy) : x(xy), y(xy) {}
	vec2_t(const T &x, const T &y) : x(x), y(y) {}
	vec2_t(const T *xy) : x(xy[0]), y(xy[1]) {}

	~vec2_t(void) {}


#pragma region Operator Overloading

	// http://en.cppreference.com/w/cpp/language/operators

#pragma region Member Access Operators

	inline T& operator[](const int offset) { return (reinterpret_cast<T*>(this))[offset]; }
	// inline T operator[](const int offset) const { return (reinterpret_cast<T*>(this))[offset]; }
	inline T operator[](const int offset) const { return ((T*) this)[offset]; }

#pragma endregion

#pragma region Arithmetic Operators

	vec2 operator+(void) const { return vec2(+this->x, +this->y); }
	vec2 operator-(void) const { return vec2(-this->x, -this->y); }

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

	vec2& operator++(void) // Prefix
	{
		++this->x;
		++this->y;

		return (*this);
	}
	inline vec2 operator++(int) { return ++(*this); } // Postfix

	vec2& operator--(void) // Prefix
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

	vec2& operator=(const vec2 &rhs)
	{
		// memcpy(this, &rhs, sizeof(rhs));
		// for (int i = 0; i < 2; i++) (*this)[i] = rhs[i];

		this->x = rhs.x;
		this->y = rhs.y;

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

	explicit inline operator T*(void) const { return reinterpret_cast<T*>(this); }

	inline operator vec2_t<float>(void) const { return vec2_t<float>(static_cast<float>(this->x), static_cast<float>(this->y)); }
	inline operator vec2_t<double>(void) const { return vec2_t<double>(static_cast<double>(this->x), static_cast<double>(this->y)); }

	inline operator vec2_t<signed int>(void) const { return vec2_t<signed int>(static_cast<signed int>(this->x), static_cast<signed int>(this->y)); }
	inline operator vec2_t<unsigned int>(void) const { return vec2_t<unsigned int>(static_cast<unsigned int>(this->x), static_cast<unsigned int>(this->y)); }

	inline operator vec2_t<bool>(void) const { return vec2_t<bool>(static_cast<bool>(this->x), static_cast<bool>(this->y)); }

	inline operator vec2_t<signed long>(void) const { return vec2_t<signed long>(static_cast<signed long>(this->x), static_cast<signed long>(this->y)); }
	inline operator vec2_t<unsigned long>(void) const { return vec2_t<unsigned long>(static_cast<unsigned long>(this->x), static_cast<unsigned long>(this->y)); }
	inline operator vec2_t<signed long long>(void) const { return vec2_t<signed long long>(static_cast<signed long long>(this->x), static_cast<signed long long>(this->y)); }
	inline operator vec2_t<unsigned long long>(void) const { return vec2_t<unsigned long long>(static_cast<unsigned long long>(this->x), static_cast<unsigned long long>(this->y)); }

#pragma endregion

#pragma region Stream Operators

#ifdef _IOSTREAM_

	friend inline std::ostream& operator<<(std::ostream &stream, const vec2 &rhs)
	{
		// return (stream << "vec2 (" << rhs.x << ", " << rhs.y << ")");
		return (stream << "vec2 {x=" << rhs.x << ", y=" << rhs.y << "}");
	}

	friend inline std::wostream& operator<<(std::wostream &stream, const vec2 &rhs)
	{
		// return (stream << L"vec2 (" << rhs.x << L", " << rhs.y << L")");
		return (stream << L"vec2 {x=" << rhs.x << L", y=" << rhs.y << L"}");
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


	// inline T magnitudeSquared(void) const { return (this->x * this->x + this->y * this->y); }
	// friend inline T magnitudeSquared(const vec2 &lhs) { return lhs.magnitudeSquared(); }

	// inline T magnitude(void) const { return sqrt(this->magnitudeSquared()); }
	// friend inline T magnitude(const vec2 &lhs) { return lhs.magnitude(); }


	inline T lengthSquared(void) const { return (this->x * this->x + this->y * this->y); }
	friend inline T lengthSquared(const vec2 &lhs) { return lhs.lengthSquared(); }

	inline T length(void) const { return sqrt(this->lengthSquared()); }
	friend inline T length(const vec2 &lhs) { return lhs.length(); }


	inline T distanceSquared(const vec2 &rhs) { return ((*this) - rhs).lengthSquared(); }
	friend inline T distanceSquared(const vec2 &lhs, const vec2 &rhs) { return lhs.distanceSquared(rhs); }

	inline T distance(const vec2 &rhs) { return ((*this) - rhs).length(); }
	friend inline T distance(const vec2 &lhs, const vec2 &rhs) { return lhs.distance(rhs); }


	vec2 normalize(const T &to = T(1.0)) const;
	friend inline vec2 normalize(const vec2 &lhs, const T &to = T(1.0)) { return lhs.normalize(to); }


	T angle(const vec2 &rhs) const
	{
		// return acos(this->dot(rhs) / (this->length() * rhs.length()));
		return (this->dot(rhs) / (this->length() * rhs.length()));
	}
	friend inline T angle(const vec2 &lhs, const vec2 &rhs) { return lhs.angle(rhs); }


	// Calculates the projection of a onto b
	// 
	// Reference: http://en.wikipedia.org/wiki/Vector_projection#Vector_projection_2
	inline vec2 project(const vec2 &b) const
	{
		const float length = b.length();

		return ((this->dot(b) / (length * length)) * b);
	}
	friend inline vec2 project(const vec2 &a, const vec2 &b) { return a.project(b); }


	// Calculates the components of a perpendicular to b
	inline vec2 perpendicular(const vec2 &b) const
	{
		const float length = b.length();

		return ((*this) - ((this->dot(b) / (length * length)) * b));
	}
	friend inline vec2 perpendicular(const vec2 &a, const vec2 &b) { return a.perpendicular(b); }


	// Calculates the reflection vector from entering ray direction a and surface normal b
	inline vec2 reflect(const vec2 &b) const
	{
		// return a - TYPE(2.0) * vec2::project(a, b);
		// return (TYPE(2.0) * vec2::project(a, b) - a);

		return (T(2) * this->project(b) - (*this));
	}
	friend inline vec2 reflect(const vec2 &a, const vec2 &b) { return a.reflect(b); }


	inline T cosine(const vec2 &b) const
	{
		return this->normalize().dot(b.normalize());
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

	inline bool isNullVector(void) const;
	friend inline bool isNullVector(const vec2 &v) { return v.isNullVector(); }

	inline bool isUnitVector(void) const;
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

	inline vec2 abs(void) const { return vec2(abs(this->x), abs(this->y), abs(this->z)); }
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
		return this->min(max).max(min);
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
		const T dot = this->normalize().dot(to.normalize());
		const T theta = acos(dot);
		const T s = sin(theta);

		return (sin((T(1) - t) * theta) / s * (*this) + sin(t * theta) / s * to);
	}
	friend inline vec2 slerp(const vec2 &from, const vec2 &to, const T &t) { return from.slerp(to, t); }

	// Reference: https://en.wikipedia.org/wiki/Slerp
	// Referemce: https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
	inline vec2 slerp(const vec2 &to, const vec2 &t) const
	{
		const T dot = this->normalize().dot(to.normalize());
		const T theta = acos(dot);
		const T s = sin(theta);

		return (sin((T(1) - t) * theta) / s * (*this) + sin(t * theta) / s * to);
	}
	friend inline vec2 slerp(const vec2 &from, const vec2 &to, const vec2 &t) { return from.slerp(to, t); }


	// Reference: https://en.wikipedia.org/wiki/Sign_function
	inline const vec2 signum(void) const
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
};


template<typename T>
class vec3_t
{
private:

	typedef vec2_t<T> vec2;
	typedef vec3_t<T> vec3;


public:

	static const vec3_t<T> zero;
	static const vec3_t<T> one;

	static const vec3_t<T> up, down;
	static const vec3_t<T> left, right;
	static const vec3_t<T> forward, backward;


public:

	// Performs Gram-Schmidt Orthogonalization on 2 basis vectors to turn them into orthonormal basis vectors
	static vec3 orthogonalize(const vec3 &a, vec3 &b)
	{
		b = b - b.project(a);
		b = b.normalize();
	}

	// Performs Gram-Schmidt Orthogonalization on 3 basis vectors to turn them into orthonormal basis vectors
	static void orthogonalize(const vec3 &a, vec3 &b, vec3 &c)
	{
		b = b - b.project(a);
		b = b.normalize();

		c = c - c.project(a) - c.project(b);
		c = c.normalize();
	}


public:

	T x, y, z;


public:

	vec3_t(void) : x(T(0)), y(T(0)), z(T(0)) {}

	// vec3_t(const vec3 &v) : x(v.x), y(v.y), z(v.z) {}

	vec3_t(const vec3_t<float> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)) {}
	vec3_t(const vec3_t<double> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)) {}
	vec3_t(const vec3_t<signed int> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)) {}
	vec3_t(const vec3_t<unsigned int> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)) {}
	vec3_t(const vec3_t<bool> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)) {}

	vec3_t(const vec3_t<signed long> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)) {}
	vec3_t(const vec3_t<unsigned long> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)) {}
	vec3_t(const vec3_t<signed long long> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)) {}
	vec3_t(const vec3_t<unsigned long long> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)) {}

	vec3_t(const T &xyz) : x(xyz), y(xyz), z(xyz) {}
	vec3_t(const T &x, const T &y, const T &z) : x(x), y(y), z(z) {}
	vec3_t(const T *xyz) : x(xyz[0]), y(xyz[1]), z(xyz[2]) {}

	vec3_t(const vec2 &xy) : x(xy.x), y(xy.y), z(T(0)) {}
	vec3_t(const vec2 &xy, const T &z) : x(xy.x), y(xy.y), z(z) {}
	vec3_t(const T &x, const vec2 &yz) : x(x), y(yz.x), z(yz.y) {}

	~vec3_t(void) {}


#pragma region Operator Overloading

	// http://en.cppreference.com/w/cpp/language/operators

#pragma region Member Access Operators

	inline T& operator[](const int offset) { return (reinterpret_cast<T*>(this))[offset]; }
	// inline T operator[](const int offset) const { return (reinterpret_cast<T*>(this))[offset]; }
	inline T operator[](const int offset) const { return ((T*) this)[offset]; }

#pragma endregion

#pragma region Arithmetic Operators

	vec3 operator+(void) const { return vec3(+this->x, +this->y, +this->z); }
	vec3 operator-(void) const { return vec3(-this->x, -this->y, -this->z); }

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

	vec3& operator++(void) // Prefix
	{
		++this->x;
		++this->y;
		++this->z;

		return (*this);
	}
	inline vec3 operator++(int) { return ++(*this); } // Postfix

	vec3& operator--(void) // Prefix
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

	vec3& operator=(const vec2 &rhs)
	{
		(*this) = vec3(rhs, this->z, this->w);

		return (*this);
	}

	vec3& operator=(const vec3 &rhs)
	{
		// memcpy(this, &rhs, sizeof(rhs));
		// for (int i = 0; i < 3; i++) (*this)[i] = rhs[i];

		this->x = rhs.x;
		this->y = rhs.y;
		this->z = rhs.z;

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

	explicit inline operator T*(void) const { return reinterpret_cast<T*>(this); }

	inline operator vec3_t<float>(void) const { return vec3_t<float>(static_cast<float>(this->x), static_cast<float>(this->y), static_cast<float>(this->z)); }
	inline operator vec3_t<double>(void) const { return vec3_t<double>(static_cast<double>(this->x), static_cast<double>(this->y), static_cast<double>(this->z)); }

	inline operator vec3_t<signed int>(void) const { return vec3_t<signed int>(static_cast<signed int>(this->x), static_cast<signed int>(this->y), static_cast<signed int>(this->z)); }
	inline operator vec3_t<unsigned int>(void) const { return vec3_t<unsigned int>(static_cast<unsigned int>(this->x), static_cast<unsigned int>(this->y), static_cast<unsigned int>(this->z)); }

	inline operator vec3_t<bool>(void) const { return vec3_t<bool>(static_cast<bool>(this->x), static_cast<bool>(this->y), static_cast<bool>(this->z)); }

	inline operator vec3_t<signed long>(void) const { return vec3_t<signed long>(static_cast<signed long>(this->x), static_cast<signed long>(this->y), static_cast<signed long>(this->z)); }
	inline operator vec3_t<unsigned long>(void) const { return vec3_t<unsigned long>(static_cast<unsigned long>(this->x), static_cast<unsigned long>(this->y), static_cast<unsigned long>(this->z)); }
	inline operator vec3_t<signed long long>(void) const { return vec3_t<signed long long>(static_cast<signed long long>(this->x), static_cast<signed long long>(this->y), static_cast<signed long long>(this->z)); }
	inline operator vec3_t<unsigned long long>(void) const { return vec3_t<unsigned long long>(static_cast<unsigned long long>(this->x), static_cast<unsigned long long>(this->y), static_cast<unsigned long long>(this->z)); }

#pragma endregion

#pragma region Stream Operators

#ifdef _IOSTREAM_

	friend inline std::ostream& operator<<(std::ostream &stream, const vec3 &rhs)
	{
		// return (stream << "vec3 (" << rhs.x << ", " << rhs.y << ", " << rhs.z << ")");
		return (stream << "vec3 {x=" << rhs.x << ", y=" << rhs.y << ", z=" << rhs.z << "}");
	}

	friend inline std::wostream& operator<<(std::wostream &stream, const vec3 &rhs)
	{
		// return (stream << L"vec3 (" << rhs.x << L", " << rhs.y << L", " << rhs.z << L")");
		return (stream << L"vec3 {x=" << rhs.x << L", y=" << rhs.y << L", z=" << rhs.z << L"}");
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


	// inline T magnitudeSquared(void) const { return (this->x * this->x + this->y * this->y + this->z * this->z); }
	// friend inline T magnitudeSquared(const vec3 &lhs) { return lhs.magnitudeSquared(); }

	// inline T magnitude(void) const { return sqrt(this->magnitudeSquared()); }
	// friend inline T magnitude(const vec3 &lhs) { return lhs.magnitude(); }


	inline T lengthSquared(void) const { return (this->x * this->x + this->y * this->y + this->z * this->z); }
	friend inline T lengthSquared(const vec3 &lhs) { return lhs.lengthSquared(); }

	inline T length(void) const { return sqrt(this->lengthSquared()); }
	friend inline T length(const vec3 &lhs) { return lhs.length(); }


	inline T distanceSquared(const vec3 &rhs) { return ((*this) - rhs).lengthSquared(); }
	friend inline T distanceSquared(const vec3 &lhs, const vec3 &rhs) { return lhs.distanceSquared(rhs); }

	inline T distance(const vec3 &rhs) { return ((*this) - rhs).length(); }
	friend inline T distance(const vec3 &lhs, const vec3 &rhs) { return lhs.distance(rhs); }


	vec3 normalize(const T &to = T(1.0)) const;
	friend inline vec3 normalize(const vec3 &lhs, const T &to = T(1.0)) { return lhs.normalize(to); }


	T angle(const vec3 &rhs) const
	{
		// return acos(this->dot(rhs) / (this->length() * rhs.length()));
		return (this->dot(rhs) / (this->length() * rhs.length()));
	}
	friend inline T angle(const vec3 &lhs, const vec3 &rhs) { return lhs.angle(rhs); }


	// Calculates the projection of a onto b
	// 
	// Reference: http://en.wikipedia.org/wiki/Vector_projection#Vector_projection_2
	inline vec3 project(const vec3 &b) const
	{
		const float length = b.length();

		return ((this->dot(b) / (length * length)) * b);
	}
	friend inline vec3 project(const vec3 &a, const vec3 &b) { return a.project(b); }


	// Calculates the components of a perpendicular to b
	inline vec3 perpendicular(const vec3 &b) const
	{
		const float length = b.length();

		return ((*this) - ((this->dot(b) / (length * length)) * b));
	}
	friend inline vec3 perpendicular(const vec3 &a, const vec3 &b) { return a.perpendicular(b); }


	// Calculates the reflection vector from entering ray direction a and surface normal b
	inline vec3 reflect(const vec3 &b) const
	{
		// return a - TYPE(2.0) * vec3::project(a, b);
		// return (TYPE(2.0) * vec3::project(a, b) - a);

		return (T(2) * this->project(b) - (*this));
	}
	friend inline vec3 reflect(const vec3 &a, const vec3 &b) { return a.reflect(b); }


	inline T cosine(const vec3 &b) const
	{
		return this->normalize().dot(b.normalize());
	}
	friend inline T cosine(const vec3 &a, const vec3 &b) { return a.cosine(b); }

#pragma endregion

#pragma region

	inline bool isNullVector(void) const;
	friend inline bool isNullVector(const vec3 &v) { return v.isNullVector(); }

	inline bool isUnitVector(void) const;
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

	inline vec3 abs(void) const { return vec3(abs(this->x), abs(this->y), abs(this->z)); }
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
		return this->min(max).max(min);
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
		const T dot = this->normalize().dot(to.normalize());
		const T theta = acos(dot);
		const T s = sin(theta);

		return (sin((T(1) - t) * theta) / s * (*this) + sin(t * theta) / s * to);
	}
	friend inline vec3 slerp(const vec3 &from, const vec3 &to, const T &t) { return from.slerp(to, t); }

	// Reference: https://en.wikipedia.org/wiki/Slerp
	// Referemce: https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
	inline vec3 slerp(const vec3 &to, const vec3 &t) const
	{
		const T dot = this->normalize().dot(to.normalize());
		const T theta = acos(dot);
		const T s = sin(theta);

		return (sin((T(1) - t) * theta) / s * (*this) + sin(t * theta) / s * to);
	}
	friend inline vec3 slerp(const vec3 &from, const vec3 &to, const vec3 &t) { return from.slerp(to, t); }


	// Reference: https://en.wikipedia.org/wiki/Sign_function
	inline const vec3 signum(void) const
	{
		return vec3(
			((this->x < T(0)) ? T(-1) : ((this->x > T(0)) ? T(1) : T(0))),
			((this->y < T(0)) ? T(-1) : ((this->y > T(0)) ? T(1) : T(0))),
			((this->z < T(0)) ? T(-1) : ((this->z > T(0)) ? T(1) : T(0)))
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

	vec4_t(void) : x(T(0)), y(T(0)), z(T(0)), w(T(0)) {}

	// vec4_t(const vec4 &v) : x(v.x), y(v.y), z(v.z), w(v.w) {}

	vec4_t(const vec4_t<float> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(T(v.w)) {}
	vec4_t(const vec4_t<double> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(T(v.w)) {}
	vec4_t(const vec4_t<signed int> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(T(v.w)) {}
	vec4_t(const vec4_t<unsigned int> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(T(v.w)) {}
	vec4_t(const vec4_t<bool> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(T(v.w)) {}

	vec4_t(const vec4_t<signed long> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(T(v.w)) {}
	vec4_t(const vec4_t<unsigned long> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(T(v.w)) {}
	vec4_t(const vec4_t<signed long long> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(T(v.w)) {}
	vec4_t(const vec4_t<unsigned long long> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(T(v.w)) {}

	vec4_t(const T &xyzw) : x(xyzw), y(xyzw), z(xyzw), w(xyzw) {}
	vec4_t(const T &x, const T &y, const T &z, const T &w) : x(x), y(y), z(z), w(w) {}
	vec4_t(const T *xyzw) : x(xyzw[0]), y(xyzw[1]), z(xyzw[2]), w(xyzw[3]) {}

	vec4_t(const vec2 &xy) : x(xy.x), y(xy.y), z(T(0)), w(T(0)) {}
	vec4_t(const vec2 &xy, const T &z, const T &w) : x(xy.x), y(xy.y), z(z), w(w) {}
	vec4_t(const T &x, const T &y, const vec2 &zw) : x(x), y(y), z(zw.x), w(zw.y) {}
	vec4_t(const vec2 &xy, const vec2 &zw) : x(xy.x), y(xy.y), z(zw.x), w(zw.y) {}

	vec4_t(const vec3 &xyz) : x(xyz.x), y(xyz.y), z(xyz.z), w(T(0)) {}
	vec4_t(const vec3 &xyz, const T &w) : x(xyz.x), y(xyz.y), z(xyz.z), w(w) {}
	vec4_t(const T &x, const vec3 &yzw) : x(x), y(yzw.x), z(yzw.y), w(yzw.z) {}

	~vec4_t(void) {}


#pragma region Operator Overloading

	// http://en.cppreference.com/w/cpp/language/operators

#pragma region Member Access Operators

	inline T& operator[](const int offset) { return (reinterpret_cast<T*>(this))[offset]; }
	// inline T operator[](const int offset) const { return (reinterpret_cast<T*>(this))[offset]; }
	inline T operator[](const int offset) const { return ((T*) this)[offset]; }

#pragma endregion

#pragma region Arithmetic Operators

	vec4 operator+(void) const { return vec4(+this->x, +this->y, +this->z, +this->w); }
	vec4 operator-(void) const { return vec4(-this->x, -this->y, -this->z, -this->w); }

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

	vec4& operator++(void) // Prefix
	{
		++this->x;
		++this->y;
		++this->z;
		++this->w;

		return (*this);
	}
	inline vec4 operator++(int) { return ++(*this); } // Postfix

	vec4& operator--(void) // Prefix
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

	vec4& operator=(const vec2 &rhs)
	{
		(*this) = vec4(rhs, this->z, this->w);

		return (*this);
	}

	vec4& operator=(const vec3 &rhs)
	{
		(*this) = vec4(rhs, this->w);

		return (*this);
	}

	vec4& operator=(const vec4 &rhs)
	{
		// memcpy(this, &rhs, sizeof(rhs));
		// for (int i = 0; i < 4; i++) (*this)[i] = rhs[i];

		this->x = rhs.x;
		this->y = rhs.y;
		this->z = rhs.z;
		this->w = rhs.w;

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

	explicit inline operator T*(void) const { return reinterpret_cast<T*>(this); }

	inline operator vec4_t<float>(void) const { return vec4_t<float>(static_cast<float>(this->x), static_cast<float>(this->y), static_cast<float>(this->z), static_cast<float>(this->w)); }
	inline operator vec4_t<double>(void) const { return vec4_t<double>(static_cast<double>(this->x), static_cast<double>(this->y), static_cast<double>(this->z), static_cast<double>(this->w)); }

	inline operator vec4_t<signed int>(void) const { return vec4_t<signed int>(static_cast<signed int>(this->x), static_cast<signed int>(this->y), static_cast<signed int>(this->z), static_cast<signed int>(this->w)); }
	inline operator vec4_t<unsigned int>(void) const { return vec4_t<unsigned int>(static_cast<unsigned int>(this->x), static_cast<unsigned int>(this->y), static_cast<unsigned int>(this->z), static_cast<unsigned int>(this->w)); }

	inline operator vec4_t<bool>(void) const { return vec4_t<bool>(static_cast<bool>(this->x), static_cast<bool>(this->y), static_cast<bool>(this->z), static_cast<bool>(this->w)); }

	inline operator vec4_t<signed long>(void) const { return vec4_t<signed long>(static_cast<signed long>(this->x), static_cast<signed long>(this->y), static_cast<signed long>(this->z), static_cast<signed long>(this->w)); }
	inline operator vec4_t<unsigned long>(void) const { return vec4_t<unsigned long>(static_cast<unsigned long>(this->x), static_cast<unsigned long>(this->y), static_cast<unsigned long>(this->z), static_cast<unsigned long>(this->w)); }
	inline operator vec4_t<signed long long>(void) const { return vec4_t<signed long long>(static_cast<signed long long>(this->x), static_cast<signed long long>(this->y), static_cast<signed long long>(this->z), static_cast<signed long long>(this->w)); }
	inline operator vec4_t<unsigned long long>(void) const { return vec4_t<unsigned long long>(static_cast<unsigned long long>(this->x), static_cast<unsigned long long>(this->y), static_cast<unsigned long long>(this->z), static_cast<unsigned long long>(this->w)); }

#pragma endregion

#pragma region Stream Operators

#ifdef _IOSTREAM_

	friend inline std::ostream& operator<<(std::ostream &stream, const vec4 &rhs)
	{
		// return (stream << "vec4 (" << rhs.x << ", " << rhs.y << ", " << rhs.z << ", " << rhs.w << ")");
		return (stream << "vec4 {x=" << rhs.x << ", y=" << rhs.y << ", z=" << rhs.z << ", w=" << rhs.w << "}");
	}

	friend inline std::wostream& operator<<(std::wostream &stream, const vec4 &rhs)
	{
		// return (stream << L"vec4 (" << rhs.x << L", " << rhs.y << L", " << rhs.z << L", " << rhs.w << L")");
		return (stream << L"vec4 {x=" << rhs.x << L", y=" << rhs.y << L", z=" << rhs.z << L", w=" << rhs.w << L"}");
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


	// A 4D vector doesn't per se have a cross product, not a feasible one at least


	// inline T magnitudeSquared(void) const { return (this->x * this->x + this->y * this->y + this->z * this->z + this->w * this->w); }
	// friend inline T magnitudeSquared(const vec4 &lhs) { return lhs.magnitudeSquared(); }

	// inline T magnitude(void) const { return sqrt(this->magnitudeSquared()); }
	// friend inline T magnitude(const vec4 &lhs) { return lhs.magnitude(); }


	inline T lengthSquared(void) const { return (this->x * this->x + this->y * this->y + this->z * this->z + this->w * this->w); }
	friend inline T lengthSquared(const vec4 &lhs) { return lhs.lengthSquared(); }

	inline T length(void) const { return sqrt(this->lengthSquared()); }
	friend inline T length(const vec4 &lhs) { return lhs.length(); }


	inline T distanceSquared(const vec4 &rhs) { return ((*this) - rhs).lengthSquared(); }
	friend inline T distanceSquared(const vec4 &lhs, const vec4 &rhs) { return lhs.distanceSquared(rhs); }

	inline T distance(const vec4 &rhs) { return ((*this) - rhs).length(); }
	friend inline T distance(const vec4 &lhs, const vec4 &rhs) { return lhs.distance(rhs); }


	vec4 normalize(const T &to = T(1.0)) const;
	friend inline vec4 normalize(const vec4 &lhs, const T &to = T(1.0)) { return lhs.normalize(to); }


	T angle(const vec4 &rhs) const
	{
		// return acos(this->dot(rhs) / (this->length() * rhs.length()));
		return (this->dot(rhs) / (this->length() * rhs.length()));
	}
	friend inline T angle(const vec4 &lhs, const vec4 &rhs) { return lhs.angle(rhs); }

#pragma endregion

#pragma region

	inline bool isNullVector(void) const;
	friend inline bool isNullVector(const vec4 &v) { return v.isNullVector(); }

	inline bool isUnitVector(void) const;
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

	inline vec4 abs(void) const { return vec4(abs(this->x), abs(this->y), abs(this->z), abs(this->w)); }
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
		return this->min(max).max(min);
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
		const T dot = this->normalize().dot(to.normalize());
		const T theta = acos(dot);
		const T s = sin(theta);

		return (sin((T(1) - t) * theta) / s * (*this) + sin(t * theta) / s * to);
	}
	friend inline vec4 slerp(const vec4 &from, const vec4 &to, const T &t) { return from.slerp(to, t); }

	// Reference: https://en.wikipedia.org/wiki/Slerp
	// Referemce: https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
	vec4 slerp(const vec4 &to, const vec4 &t) const
	{
		const T dot = this->normalize().dot(to.normalize());
		const T theta = acos(dot);
		const T s = sin(theta);

		return (sin((T(1) - t) * theta) / s * (*this) + sin(t * theta) / s * to);
	}
	friend inline vec4 slerp(const vec4 &from, const vec4 &to, const vec4 &t) { return from.slerp(to, t); }


	// Reference: https://en.wikipedia.org/wiki/Sign_function
	inline const vec4 signum(void) const
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
};


// It isn't an optimal solution, to inline all template functions that has explicit specialization.
// But it is needed if we don't want to run into "multiple definitions" compilation error.


#pragma region Static Members

template<typename T> const vec2_t<T> vec2_t<T>::zero = vec2_t<T>(T(0), T(0));
template<typename T> const vec2_t<T> vec2_t<T>::one = vec2_t<T>(T(1), T(1));

template<typename T> const vec2_t<T> vec2_t<T>::up = vec2_t<T>(T(0), T(1));
template<typename T> const vec2_t<T> vec2_t<T>::down = vec2_t<T>(T(0), T(-1));

template<typename T> const vec2_t<T> vec2_t<T>::left = vec2_t<T>(T(-1), T(0));
template<typename T> const vec2_t<T> vec2_t<T>::right = vec2_t<T>(T(1), T(0));

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
	float length = this->length();

	if (!LINALG_FEQUAL(length, 0.0f) && !LINALG_FEQUAL(length, to))
	{
		length = to / length;

		return fvec2(*this) * length;
	}

	return (*this);
}

template<> inline dvec2 dvec2::normalize(const double &to) const
{
	double length = this->length();

	if (!LINALG_DEQUAL(length, 0.0) && !LINALG_DEQUAL(length, to))
	{
		length = to / length;

		return dvec2(*this) * length;
	}

	return (*this);
}

#pragma endregion

#pragma region Check Vector Kind

template<typename T> bool vec2_t<T>::isNullVector(void) const
{
	// It's faster to just compare to 0, than to calculate the length
	return ((this->x == 0) && (this->y == 0));
}

template<> bool fvec2::isNullVector(void) const
{
	// It's faster to just compare to 0, than to calculate the length
	return (LINALG_FEQUAL(this->x, 0.0f) && LINALG_FEQUAL(this->y, 0.0f));
}

template<> bool dvec2::isNullVector(void) const
{
	// It's faster to just compare to 0, than to calculate the length
	return (LINALG_DEQUAL(this->x, 0.0) && LINALG_DEQUAL(this->y, 0.0));
}


template<> inline bool fvec2::isUnitVector() const
{
	const float length = this->length();

	return LINALG_FEQUAL(length, 1.0f);
}

template<> inline bool dvec2::isUnitVector() const
{
	const double length = this->length();

	return LINALG_DEQUAL(length, 1.0);
}


template<> bool fvec2::isNormalized(const float &to) const
{
	const float length = this->length();

	return LINALG_FEQUAL(length, to);
}

template<> bool dvec2::isNormalized(const double &to) const
{
	const double length = this->length();

	return LINALG_DEQUAL(length, to);
}


template<typename T> bool vec2_t<T>::isOrthogonalTo(const vec2 &rhs) const
{
	return (this->dot(rhs) == 0);
}

template<> bool fvec2::isOrthogonalTo(const vec2 &rhs) const
{
	const float dot = this->dot(rhs);

	return (LINALG_FEQUAL(dot, 0.0f));
}

template<> bool dvec2::isOrthogonalTo(const vec2 &rhs) const
{
	const double dot = this->dot(rhs);

	return (LINALG_DEQUAL(dot, 0.0));
}


template<typename T> bool vec2_t<T>::isPerpendicularTo(const vec2 &rhs) const
{
	return (this->dot(rhs) == 0);
}

template<> bool fvec2::isPerpendicularTo(const vec2 &rhs) const
{
	const float dot = this->dot(rhs);

	return (LINALG_FEQUAL(dot, 0.0f));
}

template<> bool dvec2::isPerpendicularTo(const vec2 &rhs) const
{
	const double dot = this->dot(rhs);

	return (LINALG_DEQUAL(dot, 0.0));
}


template<typename T> bool vec2_t<T>::isParallelTo(const vec2 &rhs) const
{
	return (this->dot(rhs) == 1);
}

template<> bool fvec2::isParallelTo(const vec2 &rhs) const
{
	const float dot = this->dot(rhs);

	return (LINALG_FEQUAL(dot, 1.0f));
}

template<> bool dvec2::isParallelTo(const vec2 &rhs) const
{
	const double dot = this->dot(rhs);

	return (LINALG_DEQUAL(dot, 1.0));
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
	float length = this->length();

	if (!LINALG_FEQUAL(length, 0.0f) && !LINALG_FEQUAL(length, to))
	{
		length = to / length;

		return fvec3(*this) * length;
	}

	return (*this);
}

template<> inline dvec3 dvec3::normalize(const double &to) const
{
	double length = this->length();

	if (!LINALG_DEQUAL(length, 0.0) && !LINALG_DEQUAL(length, to))
	{
		length = to / length;

		return dvec3(*this) * length;
	}

	return (*this);
}

#pragma endregion

#pragma region Check Vector Kind

template<typename T> bool vec3_t<T>::isNullVector(void) const
{
	// It's faster to just compare to 0, than to calculate the length
	return ((this->x == 0) && (this->y == 0) && (this->z == 0));
}

template<> bool fvec3::isNullVector(void) const
{
	// It's faster to just compare to 0, than to calculate the length
	return (LINALG_FEQUAL(this->x, 0.0f) && LINALG_FEQUAL(this->y, 0.0f) && LINALG_FEQUAL(this->z, 0.0f));
}

template<> bool dvec3::isNullVector(void) const
{
	// It's faster to just compare to 0, than to calculate the length
	return (LINALG_DEQUAL(this->x, 0.0) && LINALG_DEQUAL(this->y, 0.0) && LINALG_DEQUAL(this->z, 0.0));
}


template<> bool fvec3::isUnitVector() const
{
	const float length = this->length();

	return LINALG_FEQUAL(length, 1.0f);
}

template<> bool dvec3::isUnitVector() const
{
	const double length = this->length();

	return LINALG_DEQUAL(length, 1.0);
}


template<> bool fvec3::isNormalized(const float &to) const
{
	const float length = this->length();

	return LINALG_FEQUAL(length, to);
}

template<> bool dvec3::isNormalized(const double &to) const
{
	const double length = this->length();

	return LINALG_DEQUAL(length, to);
}


template<typename T> bool vec3_t<T>::isOrthogonalTo(const vec3 &rhs) const
{
	return (this->dot(rhs) == 0);
}

template<> bool fvec3::isOrthogonalTo(const vec3 &rhs) const
{
	const float dot = this->dot(rhs);

	return (LINALG_FEQUAL(dot, 0.0f));
}

template<> bool dvec3::isOrthogonalTo(const vec3 &rhs) const
{
	const double dot = this->dot(rhs);

	return (LINALG_DEQUAL(dot, 0.0));
}


template<typename T> bool vec3_t<T>::isPerpendicularTo(const vec3 &rhs) const
{
	return (this->dot(rhs) == 0);
}

template<> bool fvec3::isPerpendicularTo(const vec3 &rhs) const
{
	const float dot = this->dot(rhs);

	return (LINALG_FEQUAL(dot, 0.0f));
}

template<> bool dvec3::isPerpendicularTo(const vec3 &rhs) const
{
	const double dot = this->dot(rhs);

	return (LINALG_DEQUAL(dot, 0.0));
}


template<typename T> bool vec3_t<T>::isParallelTo(const vec3 &rhs) const
{
	return (this->dot(rhs) == 1);
}

template<> bool fvec3::isParallelTo(const vec3 &rhs) const
{
	const float dot = this->dot(rhs);

	return (LINALG_FEQUAL(dot, 1.0f));
}

template<> bool dvec3::isParallelTo(const vec3 &rhs) const
{
	const double dot = this->dot(rhs);

	return (LINALG_DEQUAL(dot, 1.0));
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



#pragma region Static Members

template<typename T> const vec4_t<T> vec4_t<T>::zero = vec4_t<T>(T(0), T(0), T(0), T(0));
template<typename T> const vec4_t<T> vec4_t<T>::one = vec4_t<T>(T(1), T(1), T(1), T(1));

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
	float length = this->length();

	if (!LINALG_FEQUAL(length, 0.0f) && !LINALG_FEQUAL(length, to))
	{
		length = to / length;

		return fvec4(*this) * length;
	}

	return (*this);
}

template<> inline dvec4 dvec4::normalize(const double &to) const
{
	double length = this->length();

	if (!LINALG_DEQUAL(length, 0.0) && !LINALG_DEQUAL(length, to))
	{
		length = to / length;

		return dvec4(*this) * length;
	}

	return (*this);
}

#pragma endregion

#pragma region Check Vector Kind

template<typename T> bool vec4_t<T>::isNullVector(void) const
{
	// It's faster to just compare to 0, than to calculate the length
	return ((this->x == 0) && (this->y == 0) && (this->z == 0) && (this->w == 0));
}

template<> bool fvec4::isNullVector(void) const
{
	// It's faster to just compare to 0, than to calculate the length
	return (LINALG_FEQUAL(this->x, 0.0f) && LINALG_FEQUAL(this->y, 0.0f) && LINALG_FEQUAL(this->z, 0.0f) && LINALG_FEQUAL(this->w, 0.0f));
}

template<> bool dvec4::isNullVector(void) const
{
	// It's faster to just compare to 0, than to calculate the length
	return (LINALG_DEQUAL(this->x, 0.0) && LINALG_DEQUAL(this->y, 0.0) && LINALG_DEQUAL(this->z, 0.0) && LINALG_DEQUAL(this->w, 0.0));
}


template<> bool fvec4::isUnitVector() const
{
	const float length = this->length();

	return LINALG_FEQUAL(length, 1.0f);
}

template<> bool dvec4::isUnitVector() const
{
	const double length = this->length();

	return LINALG_DEQUAL(length, 1.0);
}


template<> bool fvec4::isNormalized(const float &to) const
{
	const float length = this->length();

	return LINALG_FEQUAL(length, to);
}

template<> bool dvec4::isNormalized(const double &to) const
{
	const double length = this->length();

	return LINALG_DEQUAL(length, to);
}


template<typename T> bool vec4_t<T>::isOrthogonalTo(const vec4 &rhs) const
{
	return (this->dot(rhs) == 0);
}

template<> bool fvec4::isOrthogonalTo(const vec4 &rhs) const
{
	const float dot = this->dot(rhs);

	return (LINALG_FEQUAL(dot, 0.0f));
}

template<> bool dvec4::isOrthogonalTo(const vec4 &rhs) const
{
	const double dot = this->dot(rhs);

	return (LINALG_DEQUAL(dot, 0.0));
}


template<typename T> bool vec4_t<T>::isPerpendicularTo(const vec4 &rhs) const
{
	return (this->dot(rhs) == 0);
}

template<> bool fvec4::isPerpendicularTo(const vec4 &rhs) const
{
	const float dot = this->dot(rhs);

	return (LINALG_FEQUAL(dot, 0.0f));
}

template<> bool dvec4::isPerpendicularTo(const vec4 &rhs) const
{
	const double dot = this->dot(rhs);

	return (LINALG_DEQUAL(dot, 0.0));
}


template<typename T> bool vec4_t<T>::isParallelTo(const vec4 &rhs) const
{
	return (this->dot(rhs) == 1);
}

template<> bool fvec4::isParallelTo(const vec4 &rhs) const
{
	const float dot = this->dot(rhs);

	return (LINALG_FEQUAL(dot, 1.0f));
}

template<> bool dvec4::isParallelTo(const vec4 &rhs) const
{
	const double dot = this->dot(rhs);

	return (LINALG_DEQUAL(dot, 1.0));
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


// Enable structure padding
#pragma pack(pop)


#endif