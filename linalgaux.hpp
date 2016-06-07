
// Author: Christian Vallentin <mail@vallentinsource.com>
// Website: http://vallentinsource.com
// Repository: https://github.com/MrVallentin/LinearAlgebra
//
// Date Created: October 01, 2013
// Last Modified: June 07, 2016

#ifndef LINEAR_ALGEBRA_AUXILIARY_HPP
#define LINEAR_ALGEBRA_AUXILIARY_HPP


#include <stack>

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


	void pushMatrix()
	{
		mat4 m = mat4::identity;

		if (!stack.empty())
			m = stack.top();

		stack.push(m);
	}

	void popMatrix()
	{
		if (!stack.empty())
			stack.pop();
	}

	size_t size() const
	{
		return stack.size();
	}

	void clear()
	{
		while (!stack.empty())
			stack.pop();
	}


	void loadMatrix(const mat4 &matrix)
	{
		if (stack.empty())
			pushMatrix();

		stack.pop();
		stack.push(matrix);
	}

	void loadIdentity()
	{
		loadMatrix(mat4::identity);
	}

	void multMatrix(const mat4 &matrix)
	{
		if (stack.empty())
			pushMatrix();

		loadMatrix(stack.top() * matrix);
	}
};


#endif