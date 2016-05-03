
// Author: Christian Vallentin <mail@vallentinsource.com>
// Website: http://vallentinsource.com
// Repository: https://github.com/MrVallentin/LinearMath
//
// Date Created: January 30, 2015
// Last Modified: May 03, 2016

#ifndef LM_VEC4_H
#define LM_VEC4_H


#include <math.h>


#include "vec2.h"
#include "vec3.h"


// #ifndef FEQUAL
// #	define FEQUAL(x, y) (abs(a - b) <= 1E-6)
// #endif


#undef FEQUAL
#undef DEQUAL

// #define FEQUAL(x, y) ((((y) - 1E-6f) < (x)) && ((x) < ((y) + 1E-6f)))
// #define DEQUAL(x, y) ((((y) - 1E-6) < (x)) && ((x) < ((y) + 1E-6)))

// This was changed from 1E-6 to 1E-4 as asserting rotate(90deg) didn't equal
#define FEQUAL(x, y) ((((y) - 1E-4f) < (x)) && ((x) < ((y) + 1E-4f)))
#define DEQUAL(x, y) ((((y) - 1E-4) < (x)) && ((x) < ((y) + 1E-4)))


#ifndef vec_t
#	ifdef scalar_t
#		define vec_t scalar_t
#	else
#		define vec_t float
#	endif
#endif


// Disable structure padding
#pragma pack(push, 1)


template<typename T> class tvec4;


template<typename T> using vec4_t = tvec4<T>;

typedef tvec4<vec_t> vec4;

typedef tvec4<float> fvec4;
typedef tvec4<double> dvec4;

typedef tvec4<signed int> ivec4;
typedef tvec4<unsigned int> uvec4;

typedef tvec4<bool> bvec4;

typedef tvec4<signed long> lvec4;
typedef tvec4<unsigned long> ulvec4;

typedef tvec4<signed long long> llvec4;
typedef tvec4<unsigned long long> ullvec4;


template<typename T>
class tvec4
{
private:

	typedef tvec2<T> vec2;
	typedef tvec3<T> vec3;
	typedef tvec4<T> vec4;


public:

	static const tvec4<T> zero;
	static const tvec4<T> one;


public:

	T x, y, z, w;


public:

	tvec4(void) : x(T(0)), y(T(0)), z(T(0)), w(T(0)) {}

	// tvec4(const vec4 &v) : x(v.x), y(v.y), z(v.z), w(v.w) {}

	tvec4(const tvec4<float> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(T(v.w)) {}
	tvec4(const tvec4<double> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(T(v.w)) {}
	tvec4(const tvec4<signed int> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(T(v.w)) {}
	tvec4(const tvec4<unsigned int> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(T(v.w)) {}
	tvec4(const tvec4<bool> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(T(v.w)) {}

	tvec4(const tvec4<signed long> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(T(v.w)) {}
	tvec4(const tvec4<unsigned long> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(T(v.w)) {}
	tvec4(const tvec4<signed long long> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(T(v.w)) {}
	tvec4(const tvec4<unsigned long long> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(T(v.w)) {}

	tvec4(const T &xyzw) : x(xyzw), y(xyzw), z(xyzw), w(xyzw) {}
	tvec4(const T &x, const T &y, const T &z, const T &w) : x(x), y(y), z(z), w(w) {}
	tvec4(const T *xyzw) : x(xyzw[0]), y(xyzw[1]), z(xyzw[2]), w(xyzw[3]) {}

	tvec4(const vec2 &xy) : x(xy.x), y(xy.y), z(T(0)), w(T(0)) {}
	tvec4(const vec2 &xy, const T &z, const T &w) : x(xy.x), y(xy.y), z(z), w(w) {}
	tvec4(const T &x, const T &y, const vec2 &zw) : x(x), y(y), z(zw.x), w(zw.y) {}
	tvec4(const vec2 &xy, const vec2 &zw) : x(xy.x), y(xy.y), z(zw.x), w(zw.y) {}

	tvec4(const vec3 &xyz) : x(xyz.x), y(xyz.y), z(xyz.z), w(T(0)) {}
	tvec4(const vec3 &xyz, const T &w) : x(xyz.x), y(xyz.y), z(xyz.z), w(w) {}
	tvec4(const T &x, const vec3 &yzw) : x(x), y(yzw.x), z(yzw.y), w(yzw.z) {}

	~tvec4(void) {}


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

	// Declaration here... The definition can be found after the class.
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

	inline operator tvec4<float>(void) const { return tvec4<float>(static_cast<float>(this->x), static_cast<float>(this->y), static_cast<float>(this->z), static_cast<float>(this->w)); }
	inline operator tvec4<double>(void) const { return tvec4<double>(static_cast<double>(this->x), static_cast<double>(this->y), static_cast<double>(this->z), static_cast<double>(this->w)); }

	inline operator tvec4<signed int>(void) const { return tvec4<signed int>(static_cast<signed int>(this->x), static_cast<signed int>(this->y), static_cast<signed int>(this->z), static_cast<signed int>(this->w)); }
	inline operator tvec4<unsigned int>(void) const { return tvec4<unsigned int>(static_cast<unsigned int>(this->x), static_cast<unsigned int>(this->y), static_cast<unsigned int>(this->z), static_cast<unsigned int>(this->w)); }
	
	inline operator tvec4<bool>(void) const { return tvec4<bool>(static_cast<bool>(this->x), static_cast<bool>(this->y), static_cast<bool>(this->z), static_cast<bool>(this->w)); }

	inline operator tvec4<signed long>(void) const { return tvec4<signed long>(static_cast<signed long>(this->x), static_cast<signed long>(this->y), static_cast<signed long>(this->z), static_cast<signed long>(this->w)); }
	inline operator tvec4<unsigned long>(void) const { return tvec4<unsigned long>(static_cast<unsigned long>(this->x), static_cast<unsigned long>(this->y), static_cast<unsigned long>(this->z), static_cast<unsigned long>(this->w)); }
	inline operator tvec4<signed long long>(void) const { return tvec4<signed long long>(static_cast<signed long long>(this->x), static_cast<signed long long>(this->y), static_cast<signed long long>(this->z), static_cast<signed long long>(this->w)); }
	inline operator tvec4<unsigned long long>(void) const { return tvec4<unsigned long long>(static_cast<unsigned long long>(this->x), static_cast<unsigned long long>(this->y), static_cast<unsigned long long>(this->z), static_cast<unsigned long long>(this->w)); }

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


	// Declaration here... The definition can be found after the class.
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

	// Declaration here... The definition can be found after the class.
	inline bool isNullVector(void) const;
	friend inline bool isNullVector(const vec4 &v) { return v.isNullVector(); }

	// Declaration here... The definition can be found after the class.
	inline bool isUnitVector(void) const;
	friend inline bool isUnitVector(const vec4 &v) { return v.isUnitVector(); }

	// Declaration here... The definition can be found after the class.
	bool isNormalized(const T &to = T(1.0)) const;
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

#define SWIZZLE_INDEX(index, c) \
		if ((c == 'X') || (c == 'X') || (c == 'r') || (c == 'R') || (c == 's') || (c == 'S')) index = 0; \
		else if ((c == 'y') || (c == 'Y') || (c == 'g') || (c == 'G') || (c == 't') || (c == 'T')) index = 1; \
		else if ((c == 'z') || (c == 'Z') || (c == 'b') || (c == 'B') || (c == 'p') || (c == 'P')) index = 2; \
		else if ((c == 'w') || (c == 'W') || (c == 'a') || (c == 'A') || (c == 'q') || (c == 'Q')) index = 3; \
		else index = 0;
	
	vec2 swizzle(const char x, const char y) const
	{
		int x_index = 0, y_index = 1;

		SWIZZLE_INDEX(x_index, x);
		SWIZZLE_INDEX(y_index, y);

		return vec2((*this)[x_index], (*this)[y_index]);
	}

	vec3 swizzle(const char x, const char y, const char z) const
	{
		int x_index = 0, y_index = 1, z_index = 2;

		SWIZZLE_INDEX(x_index, x);
		SWIZZLE_INDEX(y_index, y);
		SWIZZLE_INDEX(z_index, z);

		return vec3((*this)[x_index], (*this)[y_index], (*this)[z_index]);
	}

	vec4 swizzle(const char x, const char y, const char z, const char w) const
	{
		int x_index = 0, y_index = 1, z_index = 2, w_index = 3;
		
		SWIZZLE_INDEX(x_index, x);
		SWIZZLE_INDEX(y_index, y);
		SWIZZLE_INDEX(z_index, z);
		SWIZZLE_INDEX(w_index, w);

		return vec4((*this)[x_index], (*this)[y_index], (*this)[z_index], (*this)[w_index]);
	}

#undef SWIZZLE_INDEX

#pragma endregion
};


#pragma region Static Members

template<typename T> const tvec4<T> tvec4<T>::zero = tvec4<T>(T(0), T(0), T(0), T(0));
template<typename T> const tvec4<T> tvec4<T>::one = tvec4<T>(T(1), T(1), T(1), T(1));

#pragma endregion


#pragma region Validate sizeof Templated Objects

#ifndef STATIC_ASSERT
#	define STATIC_ASSERT(bool_constexpr) static_assert(bool_constexpr, #bool_constexpr)
#endif

STATIC_ASSERT(sizeof(tvec4<float>) == (sizeof(float) * 4));
STATIC_ASSERT(sizeof(tvec4<double>) == (sizeof(double) * 4));

STATIC_ASSERT(sizeof(tvec4<signed int>) == (sizeof(signed int) * 4));
STATIC_ASSERT(sizeof(tvec4<unsigned int>) == (sizeof(unsigned int) * 4));

STATIC_ASSERT(sizeof(tvec4<bool>) == (sizeof(bool) * 4));

STATIC_ASSERT(sizeof(tvec4<signed long>) == (sizeof(signed long) * 4));
STATIC_ASSERT(sizeof(tvec4<unsigned long>) == (sizeof(unsigned long) * 4));

STATIC_ASSERT(sizeof(tvec4<signed long long>) == (sizeof(signed long long) * 4));
STATIC_ASSERT(sizeof(tvec4<unsigned long long>) == (sizeof(unsigned long long) * 4));

#pragma endregion


// Enable structure padding
#pragma pack(pop)


#endif