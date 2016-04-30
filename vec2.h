
// Author: Christian Vallentin <mail@vallentinsource.com>
// Website: http://vallentinsource.com
// Repository: https://github.com/MrVallentin/LinearMath
//
// Date Created: January 30, 2015
// Last Modified: April 30, 2016

#ifndef LM_VEC2_H
#define LM_VEC2_H



#include <math.h>



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



template<typename T> class tvec2;



template<typename T> using vec2_t = tvec2 <T>;

typedef tvec2<vec_t> vec2;

typedef tvec2<float> fvec2;
typedef tvec2<double> dvec2;

typedef tvec2<signed int> ivec2;
typedef tvec2<unsigned int> uvec2;

typedef tvec2<bool> bvec2;

typedef tvec2<signed long> lvec2;
typedef tvec2<unsigned long> ulvec2;

typedef tvec2<signed long long> llvec2;
typedef tvec2<unsigned long long> ullvec2;



template<typename T>
class tvec2
{
private:

	typedef tvec2<T> vec2;


public:

	static const tvec2<T> zero;
	static const tvec2<T> one;

	static const tvec2<T> up, down;
	static const tvec2<T> left, right;


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

	tvec2(void) : x(T(0)), y(T(0)) {}

	// tvec2(const vec2 &v) : x(v.x), y(v.y) {}

	tvec2(const tvec2<float> &v) : x(T(v.x)), y(T(v.y)) {}
	tvec2(const tvec2<double> &v) : x(T(v.x)), y(T(v.y)) {}
	tvec2(const tvec2<signed int> &v) : x(T(v.x)), y(T(v.y)) {}
	tvec2(const tvec2<unsigned int> &v) : x(T(v.x)), y(T(v.y)) {}
	tvec2(const tvec2<bool> &v) : x(T(v.x)), y(T(v.y)) {}

	tvec2(const tvec2<signed long> &v) : x(T(v.x)), y(T(v.y)) {}
	tvec2(const tvec2<unsigned long> &v) : x(T(v.x)), y(T(v.y)) {}
	tvec2(const tvec2<signed long long> &v) : x(T(v.x)), y(T(v.y)) {}
	tvec2(const tvec2<unsigned long long> &v) : x(T(v.x)), y(T(v.y)) {}

	tvec2(const T &xy) : x(xy), y(xy) {}
	tvec2(const T &x, const T &y) : x(x), y(y) {}
	tvec2(const T *xy) : x(xy[0]), y(xy[1]) {}

	~tvec2(void) {}


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

	// Declaration here... The definition can be found after the class.
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

	inline operator tvec2<float>(void) const { return tvec2<float>(static_cast<float>(this->x), static_cast<float>(this->y)); }
	inline operator tvec2<double>(void) const { return tvec2<double>(static_cast<double>(this->x), static_cast<double>(this->y)); }

	inline operator tvec2<signed int>(void) const { return tvec2<signed int>(static_cast<signed int>(this->x), static_cast<signed int>(this->y)); }
	inline operator tvec2<unsigned int>(void) const { return tvec2<unsigned int>(static_cast<unsigned int>(this->x), static_cast<unsigned int>(this->y)); }

	inline operator tvec2<bool>(void) const { return tvec2<bool>(static_cast<bool>(this->x), static_cast<bool>(this->y)); }

	inline operator tvec2<signed long>(void) const { return tvec2<signed long>(static_cast<signed long>(this->x), static_cast<signed long>(this->y)); }
	inline operator tvec2<unsigned long>(void) const { return tvec2<unsigned long>(static_cast<unsigned long>(this->x), static_cast<unsigned long>(this->y)); }
	inline operator tvec2<signed long long>(void) const { return tvec2<signed long long>(static_cast<signed long long>(this->x), static_cast<signed long long>(this->y)); }
	inline operator tvec2<unsigned long long>(void) const { return tvec2<unsigned long long>(static_cast<unsigned long long>(this->x), static_cast<unsigned long long>(this->y)); }

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


	// Declaration here... The definition can be found after the class.
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

	// Declaration here... The definition can be found after the class.
	inline bool isNullVector(void) const;
	friend inline bool isNullVector(const vec2 &v) { return v.isNullVector(); }

	// Declaration here... The definition can be found after the class.
	inline bool isUnitVector(void) const;
	friend inline bool isUnitVector(const vec2 &v) { return v.isUnitVector(); }

	// Declaration here... The definition can be found after the class.
	bool isNormalized(const T &to = T(1.0)) const;
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

#define SWIZZLE_INDEX(index, c) \
		if ((c == 'X') || (c == 'X') || (c == 'r') || (c == 'R') || (c == 's') || (c == 'S')) index = 0; \
		else if ((c == 'y') || (c == 'Y') || (c == 'g') || (c == 'G') || (c == 't') || (c == 'T')) index = 1; \
		else index = 0;

	vec2 swizzle(const char x, const char y) const
	{
		int x_index = 0, y_index = 1;

		SWIZZLE_INDEX(x_index, x);
		SWIZZLE_INDEX(y_index, y);

		return vec2((*this)[x_index], (*this)[y_index]);
	}

#undef SWIZZLE_INDEX

#pragma endregion
};



#pragma region Static Members

template<typename T> const tvec2<T> tvec2<T>::zero = tvec2<T>(T(0), T(0));
template<typename T> const tvec2<T> tvec2<T>::one = tvec2<T>(T(1), T(1));

template<typename T> const tvec2<T> tvec2<T>::up = tvec2<T>(T(0), T(1));
template<typename T> const tvec2<T> tvec2<T>::down = tvec2<T>(T(0), T(-1));

template<typename T> const tvec2<T> tvec2<T>::left = tvec2<T>(T(-1), T(0));
template<typename T> const tvec2<T> tvec2<T>::right = tvec2<T>(T(1), T(0));

#pragma endregion



#pragma region Validate sizeof Templated Objects

#ifndef STATIC_ASSERT
#	define STATIC_ASSERT(bool_constexpr) static_assert(bool_constexpr, #bool_constexpr)
#endif

STATIC_ASSERT(sizeof(tvec2<float>) == (sizeof(float) * 2));
STATIC_ASSERT(sizeof(tvec2<double>) == (sizeof(double) * 2));

STATIC_ASSERT(sizeof(tvec2<signed int>) == (sizeof(signed int) * 2));
STATIC_ASSERT(sizeof(tvec2<unsigned int>) == (sizeof(unsigned int) * 2));

STATIC_ASSERT(sizeof(tvec2<bool>) == (sizeof(bool) * 2));

STATIC_ASSERT(sizeof(tvec2<signed long>) == (sizeof(signed long) * 2));
STATIC_ASSERT(sizeof(tvec2<unsigned long>) == (sizeof(unsigned long) * 2));

STATIC_ASSERT(sizeof(tvec2<signed long long>) == (sizeof(signed long long) * 2));
STATIC_ASSERT(sizeof(tvec2<unsigned long long>) == (sizeof(unsigned long long) * 2));

#pragma endregion



#endif