
// Author: Christian Vallentin <mail@vallentinsource.com>
// Website: http://vallentinsource.com
// Repository: https://github.com/MrVallentin/LinearMath
//
// Date Created: January 30, 2015
// Last Modified: April 30, 2016

#include "vec4.h"



#pragma region Comparison Operators

template<typename T> bool tvec4<T>::operator==(const vec4 &rhs) const
{
	return ((this->x == rhs.x) && (this->y == rhs.y) && (this->z == rhs.z) && (this->w == rhs.w));
}

template<> bool fvec4::operator==(const fvec4 &rhs) const
{
	return (FEQUAL(this->x, rhs.x) && FEQUAL(this->y, rhs.y) && FEQUAL(this->z, rhs.z) && FEQUAL(this->w, rhs.w));
}

template<> bool dvec4::operator==(const dvec4 &rhs) const
{
	return (DEQUAL(this->x, rhs.x) && DEQUAL(this->y, rhs.y) && DEQUAL(this->z, rhs.z) && DEQUAL(this->w, rhs.w));
}

#pragma endregion

#pragma region Normalize

template<> fvec4 fvec4::normalize(const float &to) const
{
	float length = this->length();

	if (!FEQUAL(length, 0.0f) && !FEQUAL(length, to))
	{
		length = to / length;

		return fvec4(*this) * length;
	}

	return (*this);
}

template<> dvec4 dvec4::normalize(const double &to) const
{
	double length = this->length();

	if (!DEQUAL(length, 0.0) && !DEQUAL(length, to))
	{
		length = to / length;

		return dvec4(*this) * length;
	}

	return (*this);
}

#pragma endregion



#pragma region

template<typename T> bool tvec4<T>::isNullVector(void) const
{
	// It's faster to just compare to 0, than to calculate the length
	return ((this->x == 0) && (this->y == 0) && (this->z == 0) && (this->w == 0));
}

template<> bool fvec4::isNullVector(void) const
{
	// It's faster to just compare to 0, than to calculate the length
	return (FEQUAL(this->x, 0.0f) && FEQUAL(this->y, 0.0f) && FEQUAL(this->z, 0.0f) && FEQUAL(this->w, 0.0f));
}

template<> bool dvec4::isNullVector(void) const
{
	// It's faster to just compare to 0, than to calculate the length
	return (DEQUAL(this->x, 0.0) && DEQUAL(this->y, 0.0) && DEQUAL(this->z, 0.0) && DEQUAL(this->w, 0.0));
}


template<> bool fvec4::isUnitVector() const
{
	const float length = this->length();

	return FEQUAL(length, 1.0f);
}

template<> bool dvec4::isUnitVector() const
{
	const double length = this->length();

	return DEQUAL(length, 1.0);
}


template<> bool fvec4::isNormalized(const float &to) const
{
	const float length = this->length();

	return FEQUAL(length, to);
}

template<> bool dvec4::isNormalized(const double &to) const
{
	const double length = this->length();

	return DEQUAL(length, to);
}


template<typename T> bool tvec4<T>::isOrthogonalTo(const vec4 &rhs) const
{
	return (this->dot(rhs) == 0);
}

template<> bool fvec4::isOrthogonalTo(const vec4 &rhs) const
{
	const float dot = this->dot(rhs);

	return (FEQUAL(dot, 0.0f));
}

template<> bool dvec4::isOrthogonalTo(const vec4 &rhs) const
{
	const double dot = this->dot(rhs);

	return (DEQUAL(dot, 0.0));
}


template<typename T> bool tvec4<T>::isPerpendicularTo(const vec4 &rhs) const
{
	return (this->dot(rhs) == 0);
}

template<> bool fvec4::isPerpendicularTo(const vec4 &rhs) const
{
	const float dot = this->dot(rhs);

	return (FEQUAL(dot, 0.0f));
}

template<> bool dvec4::isPerpendicularTo(const vec4 &rhs) const
{
	const double dot = this->dot(rhs);

	return (DEQUAL(dot, 0.0));
}


template<typename T> bool tvec4<T>::isParallelTo(const vec4 &rhs) const
{
	return (this->dot(rhs) == 1);
}

template<> bool fvec4::isParallelTo(const vec4 &rhs) const
{
	const float dot = this->dot(rhs);

	return (FEQUAL(dot, 1.0f));
}

template<> bool dvec4::isParallelTo(const vec4 &rhs) const
{
	const double dot = this->dot(rhs);

	return (DEQUAL(dot, 1.0));
}

#pragma endregion


