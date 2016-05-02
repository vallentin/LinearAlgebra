
// Author: Christian Vallentin <mail@vallentinsource.com>
// Website: http://vallentinsource.com
// Repository: https://github.com/MrVallentin/LinearMath
//
// Date Created: January 30, 2015
// Last Modified: April 30, 2016

#include "vec2.h"



#pragma region Comparison Operators

template<typename T> bool tvec2<T>::operator==(const vec2 &rhs) const
{
	return ((this->x == rhs.x) && (this->y == rhs.y));
}

template<> bool fvec2::operator==(const fvec2 &rhs) const
{
	return (FEQUAL(this->x, rhs.x) && FEQUAL(this->y, rhs.y));
}

template<> bool dvec2::operator==(const dvec2 &rhs) const
{
	return (DEQUAL(this->x, rhs.x) && DEQUAL(this->y, rhs.y));
}

#pragma endregion

#pragma region Normalize

template<> fvec2 fvec2::normalize(const float &to) const
{
	float length = this->length();

	if (!FEQUAL(length, 0.0f) && !FEQUAL(length, to))
	{
		length = to / length;

		return fvec2(*this) * length;
	}

	return (*this);
}

template<> dvec2 dvec2::normalize(const double &to) const
{
	double length = this->length();

	if (!DEQUAL(length, 0.0) && !DEQUAL(length, to))
	{
		length = to / length;

		return dvec2(*this) * length;
	}

	return (*this);
}

#pragma endregion



#pragma region

template<typename T> bool tvec2<T>::isNullVector(void) const
{
	// It's faster to just compare to 0, than to calculate the length
	return ((this->x == 0) && (this->y == 0));
}

template<> bool fvec2::isNullVector(void) const
{
	// It's faster to just compare to 0, than to calculate the length
	return (FEQUAL(this->x, 0.0f) && FEQUAL(this->y, 0.0f));
}

template<> bool dvec2::isNullVector(void) const
{
	// It's faster to just compare to 0, than to calculate the length
	return (DEQUAL(this->x, 0.0) && DEQUAL(this->y, 0.0));
}


template<> bool fvec2::isUnitVector() const
{
	const float length = this->length();

	return FEQUAL(length, 1.0f);
}

template<> bool dvec2::isUnitVector() const
{
	const double length = this->length();

	return DEQUAL(length, 1.0);
}


template<> bool fvec2::isNormalized(const float &to) const
{
	const float length = this->length();

	return FEQUAL(length, to);
}

template<> bool dvec2::isNormalized(const double &to) const
{
	const double length = this->length();

	return DEQUAL(length, to);
}


template<typename T> bool tvec2<T>::isOrthogonalTo(const vec2 &rhs) const
{
	return (this->dot(rhs) == 0);
}

template<> bool fvec2::isOrthogonalTo(const vec2 &rhs) const
{
	const float dot = this->dot(rhs);

	return (FEQUAL(dot, 0.0f));
}

template<> bool dvec2::isOrthogonalTo(const vec2 &rhs) const
{
	const double dot = this->dot(rhs);

	return (DEQUAL(dot, 0.0));
}


template<typename T> bool tvec2<T>::isPerpendicularTo(const vec2 &rhs) const
{
	return (this->dot(rhs) == 0);
}

template<> bool fvec2::isPerpendicularTo(const vec2 &rhs) const
{
	const float dot = this->dot(rhs);

	return (FEQUAL(dot, 0.0f));
}

template<> bool dvec2::isPerpendicularTo(const vec2 &rhs) const
{
	const double dot = this->dot(rhs);

	return (DEQUAL(dot, 0.0));
}


template<typename T> bool tvec2<T>::isParallelTo(const vec2 &rhs) const
{
	return (this->dot(rhs) == 1);
}

template<> bool fvec2::isParallelTo(const vec2 &rhs) const
{
	const float dot = this->dot(rhs);

	return (FEQUAL(dot, 1.0f));
}

template<> bool dvec2::isParallelTo(const vec2 &rhs) const
{
	const double dot = this->dot(rhs);

	return (DEQUAL(dot, 1.0));
}

#pragma endregion


