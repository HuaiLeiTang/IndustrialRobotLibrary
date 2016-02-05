#include "Function.h"

#include <iostream>

namespace rovin
{
	VectorX AffineFunction::func(const VectorX & x) const
	{
		return A*x + b;
	}

	MatrixX AffineFunction::Jacobian(const VectorX & x) const
	{
		return A;
	}

	std::vector<MatrixX> AffineFunction::Hessian(const VectorX & x) const
	{
		std::vector<MatrixX> Hess(b.size());
		for (int i = 0; i < b.size(); i++)
			Hess[i] = MatrixX::Zero(x.size(), x.size());
		return Hess;
	}

	FunctionPtr AffineFunction::MultiplyConst(const Real& w) const
	{
		FunctionPtr p_wf;
		p_wf = FunctionPtr(new AffineFunction);
		std::static_pointer_cast<AffineFunction>(p_wf)->A = A*w;
		std::static_pointer_cast<AffineFunction>(p_wf)->b = b*w;
		return p_wf;
	}
	
	VectorX QuadraticFunction::func(const VectorX & x) const
	{
		return x.transpose()*A*x + b.transpose()*x + c;
	}

	MatrixX QuadraticFunction::Jacobian(const VectorX & x) const
	{
		return (A+A.transpose())*x + b;
	}

	std::vector<MatrixX> QuadraticFunction::Hessian(const VectorX & x) const
	{
		std::vector<MatrixX> Hess(1);
		Hess[0] = A + A.transpose();
		return Hess;
	}

	FunctionPtr QuadraticFunction::MultiplyConst(const Real& w) const
	{
		FunctionPtr p_wf;
		p_wf = FunctionPtr(new QuadraticFunction);
		std::static_pointer_cast<QuadraticFunction>(p_wf)->A = A*w;
		std::static_pointer_cast<QuadraticFunction>(p_wf)->b = b*w;
		std::static_pointer_cast<QuadraticFunction>(p_wf)->c = c*w;
		return p_wf;
	}

}