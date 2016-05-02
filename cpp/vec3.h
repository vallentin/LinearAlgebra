
// Author: Christian Vallentin <mail@vallentinsource.com>
// Website: http://vallentinsource.com
// Repository: https://github.com/MrVallentin/LinearMath
//
// Date Created: January 30, 2015
// Last Modified: April 30, 2016

#ifndef LM_VEC3_H
#define LM_VEC3_H



#include <math.h>



#include "vec2.h"



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



template<typename T> class tvec3;



template<typename T> using vec3_t = tvec3<T>;

typedef tvec3<vec_t> vec3;

typedef tvec3<float> fvec3;
typedef tvec3<double> dvec3;

typedef tvec3<signed int> ivec3;
typedef tvec3<unsigned int> uvec3;

typedef tvec3<bool> bvec3;

typedef tvec3<signed long> lvec3;
typedef tvec3<unsigned long> ulvec3;

typedef tvec3<signed long long> llvec3;
typedef tvec3<unsigned long long> ullvec3;



template<typename T>
class tvec3
{
private:

	typedef tvec2<T> vec2;
	typedef tvec3<T> vec3;


public:

	static const tvec3<T> zero;
	static const tvec3<T> one;

	static const tvec3<T> up, down;
	static const tvec3<T> left, right;
	static const tvec3<T> forward, backward;


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

	tvec3(void) : x(T(0)), y(T(0)), z(T(0)) {}

	// tvec3(const vec3 &v) : x(v.x), y(v.y), z(v.z) {}

	tvec3(const tvec3<float> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)) {}
	tvec3(const tvec3<double> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)) {}
	tvec3(const tvec3<signed int> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)) {}
	tvec3(const tvec3<unsigned int> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)) {}
	tvec3(const tvec3<bool> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)) {}

	tvec3(const tvec3<signed long> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)) {}
	tvec3(const tvec3<unsigned long> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)) {}
	tvec3(const tvec3<signed long long> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)) {}
	tvec3(const tvec3<unsigned long long> &v) : x(T(v.x)), y(T(v.y)), z(T(v.z)) {}

	tvec3(const T &xyz) : x(xyz), y(xyz), z(xyz) {}
	tvec3(const T &x, const T &y, const T &z) : x(x), y(y), z(z) {}
	tvec3(const T *xyz) : x(xyz[0]), y(xyz[1]), z(xyz[2]) {}

	tvec3(const vec2 &xy) : x(xy.x), y(xy.y), z(T(0)) {}
	tvec3(const vec2 &xy, const T &z) : x(xy.x), y(xy.y), z(z) {}
	tvec3(const T &x, const vec2 &yz) : x(x), y(yz.x), z(yz.y) {}

	~tvec3(void) {}


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

	// Declaration here... The definition can be found after the class.
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

	inline operator tvec3<float>(void) const { return tvec3<float>(static_cast<float>(this->x), static_cast<float>(this->y), static_cast<float>(this->z)); }
	inline operator tvec3<double>(void) const { return tvec3<double>(static_cast<double>(this->x), static_cast<double>(this->y), static_cast<double>(this->z)); }

	inline operator tvec3<signed int>(void) const { return tvec3<signed int>(static_cast<signed int>(this->x), static_cast<signed int>(this->y), static_cast<signed int>(this->z)); }
	inline operator tvec3<unsigned int>(void) const { return tvec3<unsigned int>(static_cast<unsigned int>(this->x), static_cast<unsigned int>(this->y), static_cast<unsigned int>(this->z)); }

	inline operator tvec3<bool>(void) const { return tvec3<bool>(static_cast<bool>(this->x), static_cast<bool>(this->y), static_cast<bool>(this->z)); }

	inline operator tvec3<signed long>(void) const { return tvec3<signed long>(static_cast<signed long>(this->x), static_cast<signed long>(this->y), static_cast<signed long>(this->z)); }
	inline operator tvec3<unsigned long>(void) const { return tvec3<unsigned long>(static_cast<unsigned long>(this->x), static_cast<unsigned long>(this->y), static_cast<unsigned long>(this->z)); }
	inline operator tvec3<signed long long>(void) const { return tvec3<signed long long>(static_cast<signed long long>(this->x), static_cast<signed long long>(this->y), static_cast<signed long long>(this->z)); }
	inline operator tvec3<unsigned long long>(void) const { return tvec3<unsigned long long>(static_cast<unsigned long long>(this->x), static_cast<unsigned long long>(this->y), static_cast<unsigned long long>(this->z)); }

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


	// Declaration here... The definition can be found after the class.
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

	// Declaration here... The definition can be found after the class.
	inline bool isNullVector(void) const;
	friend inline bool isNullVector(const vec3 &v) { return v.isNullVector(); }

	// Declaration here... The definition can be found after the class.
	inline bool isUnitVector(void) const;
	friend inline bool isUnitVector(const vec3 &v) { return v.isUnitVector(); }

	// Declaration here... The definition can be found after the class.
	bool isNormalized(const T &to = T(1.0)) const;
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

#define SWIZZLE_INDEX(index, c) \
		if ((c == 'X') || (c == 'X') || (c == 'r') || (c == 'R') || (c == 's') || (c == 'S')) index = 0; \
		else if ((c == 'y') || (c == 'Y') || (c == 'g') || (c == 'G') || (c == 't') || (c == 'T')) index = 1; \
		else if ((c == 'z') || (c == 'Z') || (c == 'b') || (c == 'B') || (c == 'p') || (c == 'P')) index = 2; \
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

#undef SWIZZLE_INDEX

#pragma endregion
};



#pragma region Static Members

template<typename T> const tvec3<T> tvec3<T>::zero = tvec3<T>(T(0), T(0), T(0));
template<typename T> const tvec3<T> tvec3<T>::one = tvec3<T>(T(1), T(1), T(1));

template<typename T> const tvec3<T> tvec3<T>::up = tvec3<T>(T(0), T(1), T(0));
template<typename T> const tvec3<T> tvec3<T>::down = tvec3<T>(T(0), T(-1), T(0));

template<typename T> const tvec3<T> tvec3<T>::left = tvec3<T>(T(-1), T(0), T(0));
template<typename T> const tvec3<T> tvec3<T>::right = tvec3<T>(T(1), T(0), T(0));

template<typename T> const tvec3<T> tvec3<T>::forward = tvec3<T>(T(0), T(0), T(1));
template<typename T> const tvec3<T> tvec3<T>::backward = tvec3<T>(T(0), T(0), T(-1));

#pragma endregion



#pragma region Validate sizeof Templated Objects

#ifndef STATIC_ASSERT
#	define STATIC_ASSERT(bool_constexpr) static_assert(bool_constexpr, #bool_constexpr)
#endif

STATIC_ASSERT(sizeof(tvec3<float>) == (sizeof(float) * 3));
STATIC_ASSERT(sizeof(tvec3<double>) == (sizeof(double) * 3));

STATIC_ASSERT(sizeof(tvec3<signed int>) == (sizeof(signed int) * 3));
STATIC_ASSERT(sizeof(tvec3<unsigned int>) == (sizeof(unsigned int) * 3));

STATIC_ASSERT(sizeof(tvec3<bool>) == (sizeof(bool) * 3));

STATIC_ASSERT(sizeof(tvec3<signed long>) == (sizeof(signed long) * 3));
STATIC_ASSERT(sizeof(tvec3<unsigned long>) == (sizeof(unsigned long) * 3));

STATIC_ASSERT(sizeof(tvec3<signed long long>) == (sizeof(signed long long) * 3));
STATIC_ASSERT(sizeof(tvec3<unsigned long long>) == (sizeof(unsigned long long) * 3));

#pragma endregion



#endif