
// Author: Christian Vallentin <mail@vallentinsource.com>
// Website: http://vallentinsource.com
// Repository: https://github.com/MrVallentin/LinearMath
//
// Date Created: January 30, 2015
// Last Modified: April 30, 2016

#include "vec3.h"



#pragma region Comparison Operators

template<typename T> bool tvec3<T>::operator==(const vec3 &rhs) const
{
	return ((this->x == rhs.x) && (this->y == rhs.y) && (this->z == rhs.z));
}

template<> bool fvec3::operator==(const fvec3 &rhs) const
{
	return (FEQUAL(this->x, rhs.x) && FEQUAL(this->y, rhs.y) && FEQUAL(this->z, rhs.z));
}

template<> bool dvec3::operator==(const dvec3 &rhs) const
{
	return (DEQUAL(this->x, rhs.x) && DEQUAL(this->y, rhs.y) && DEQUAL(this->z, rhs.z));
}

#pragma endregion

#pragma region Normalize

template<> fvec3 fvec3::normalize(const float &to) const
{
	float length = this->length();

	if (!FEQUAL(length, 0.0f) && !FEQUAL(length, to))
	{
		length = to / length;

		return fvec3(*this) * length;
	}

	return (*this);
}

template<> dvec3 dvec3::normalize(const double &to) const
{
	double length = this->length();

	if (!DEQUAL(length, 0.0) && !DEQUAL(length, to))
	{
		length = to / length;

		return dvec3(*this) * length;
	}

	return (*this);
}

#pragma endregion



#pragma region

template<typename T> bool tvec3<T>::isNullVector(void) const
{
	// It's faster to just compare to 0, than to calculate the length
	return ((this->x == 0) && (this->y == 0) && (this->z == 0));
}

template<> bool fvec3::isNullVector(void) const
{
	// It's faster to just compare to 0, than to calculate the length
	return (FEQUAL(this->x, 0.0f) && FEQUAL(this->y, 0.0f) && FEQUAL(this->z, 0.0f));
}

template<> bool dvec3::isNullVector(void) const
{
	// It's faster to just compare to 0, than to calculate the length
	return (DEQUAL(this->x, 0.0) && DEQUAL(this->y, 0.0) && DEQUAL(this->z, 0.0));
}


template<> bool fvec3::isUnitVector() const
{
	const float length = this->length();

	return FEQUAL(length, 1.0f);
}

template<> bool dvec3::isUnitVector() const
{
	const double length = this->length();

	return DEQUAL(length, 1.0);
}


template<> bool fvec3::isNormalized(const float &to) const
{
	const float length = this->length();

	return FEQUAL(length, to);
}

template<> bool dvec3::isNormalized(const double &to) const
{
	const double length = this->length();

	return DEQUAL(length, to);
}


template<typename T> bool tvec3<T>::isOrthogonalTo(const vec3 &rhs) const
{
	return (this->dot(rhs) == 0);
}

template<> bool fvec3::isOrthogonalTo(const vec3 &rhs) const
{
	const float dot = this->dot(rhs);

	return (FEQUAL(dot, 0.0f));
}

template<> bool dvec3::isOrthogonalTo(const vec3 &rhs) const
{
	const double dot = this->dot(rhs);

	return (DEQUAL(dot, 0.0));
}


template<typename T> bool tvec3<T>::isPerpendicularTo(const vec3 &rhs) const
{
	return (this->dot(rhs) == 0);
}

template<> bool fvec3::isPerpendicularTo(const vec3 &rhs) const
{
	const float dot = this->dot(rhs);

	return (FEQUAL(dot, 0.0f));
}

template<> bool dvec3::isPerpendicularTo(const vec3 &rhs) const
{
	const double dot = this->dot(rhs);

	return (DEQUAL(dot, 0.0));
}


template<typename T> bool tvec3<T>::isParallelTo(const vec3 &rhs) const
{
	return (this->dot(rhs) == 1);
}

template<> bool fvec3::isParallelTo(const vec3 &rhs) const
{
	const float dot = this->dot(rhs);

	return (FEQUAL(dot, 1.0f));
}

template<> bool dvec3::isParallelTo(const vec3 &rhs) const
{
	const double dot = this->dot(rhs);

	return (DEQUAL(dot, 1.0));
}

#pragma endregion


