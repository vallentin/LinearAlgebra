
// Author: Christian Vallentin <mail@vallentinsource.com>
// Website: http://vallentinsource.com
// Repository: https://github.com/MrVallentin/LinearAlgebra
// License: https://github.com/MrVallentin/LinearAlgebra/blob/master/LICENSE
//
// Date Created: October 01, 2013
// Last Modified: June 27, 2016

// Refrain from using any exposed macros, functions
// or structs prefixed with an underscore. As these
// are only intended for internal purposes. Which
// additionally means they can be removed, renamed
// or changed between minor updates without notice.

#ifndef LINEAR_ALGEBRA_AUXILIARY_HPP
#define LINEAR_ALGEBRA_AUXILIARY_HPP


#include <stack>

#include <cstddef>

#include "linalg.hpp"


template<typename T> class MatrixStackT;

typedef MatrixStackT<float> MatrixStack;
typedef MatrixStackT<double> MatrixStackD;


template<typename T>
class MatrixStackT
{
private:

	typedef mat4_t<T> mat4;


private:

	std::stack<mat4> stack;


public:

	MatrixStackT()
	{

	}

	~MatrixStackT()
	{
		clear();
	}


	mat4 getMatrix() const
	{
		if (this->stack.empty())
			return mat4::identity;

		return this->stack.top();
	}


	void pushMatrix()
	{
		this->stack.push(getMatrix());
	}

	void popMatrix()
	{
		if (!this->stack.empty())
			this->stack.pop();
	}

	size_t size() const
	{
		return this->stack.size();
	}

	void clear()
	{
		while (!this->stack.empty())
			this->stack.pop();
	}


	void loadMatrix(const mat4 &matrix)
	{
		if (this->stack.empty())
			pushMatrix();

		this->stack.pop();
		this->stack.push(matrix);
	}

	void loadIdentity()
	{
		loadMatrix(mat4::identity);
	}

	void multMatrix(const mat4 &matrix)
	{
		if (this->stack.empty())
			pushMatrix();

		loadMatrix(this->stack.top() * matrix);
	}
};


#endif