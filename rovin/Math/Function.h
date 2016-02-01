/*!
*	\file	Function.h
*	\date	2016.01.28
*	\author	Wooyoung Kim(wykim1989@gmail.com)
*	\brief	Functions
*/

#pragma once

#include "Constant.h"
#include <memory>
#include <vector>

namespace rovin
{
	class Function;
	typedef std::shared_ptr< Function > FunctionPtr;

	class Function
	{
	public:
		Function(const Real& eps = (1e-5))
		{
			_eps = eps;
		}

		VectorX operator()(const VectorX& x) const
		{
			return func(x);
		}

		virtual VectorX func(const VectorX& x) const = 0;
		virtual MatrixX Jacobian(const VectorX& x) const = 0;
		virtual std::vector< MatrixX > Hessian(const VectorX& x) const = 0;
		virtual FunctionPtr MultiplyConst(const Real& w) const = 0;
	private:
		Real _eps;
	};


	class AffineFunction : public Function
	{

	public:
		AffineFunction() {}
		MatrixX A;
		VectorX b;

		VectorX func(const VectorX& x) const;
		MatrixX Jacobian(const VectorX& x) const;
		std::vector< MatrixX > Hessian(const VectorX& x) const;
		FunctionPtr MultiplyConst(const Real& w) const;

	};

	class QuadraticFunction : public Function
	{
	public:
		QuadraticFunction() {}
		MatrixX A;
		VectorX b;
		VectorX c;

		VectorX func(const VectorX& x) const;
		MatrixX Jacobian(const VectorX& x) const;
		std::vector< MatrixX > Hessian(const VectorX& x) const;
		FunctionPtr MultiplyConst(const Real& w) const;

	};


}