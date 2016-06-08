#pragma once

#include "rovin\Math\Constant.h"
#include "rovin\Math\Common.h"
#include <iostream>


using namespace rovin;

namespace irLib
{
	namespace Opt
	{
		static void LinearConjugateGradient(const MatrixX& A, const VectorX& b, const VectorX& x0, VectorX& x)
		{
			int n = A.rows();
			int maxIter = n * 2;
			int k = 0;
			bool iterSwi = false;

			VectorX r(n), r1(n), p(n), Ap(n);
			Real alpha = 0, beta = 0;
			Real r1inner = 0;
			
			// initialization
			r = A * x0 - b;
			p = -r; x = x0;

			while (k <= maxIter)
			{
				Ap = A * p;
				alpha = (r.transpose() * r);
				alpha /= (p.transpose()*Ap);
				x += alpha*p;
				r1 = r + alpha * Ap;
				r1inner = (r1.transpose()*r1);
				beta = r1inner / (r.transpose()*r);
				p = -r1 + beta*p;

				r = r1;
				if(r1inner < 1E-10)
				{
					iterSwi = true;
					break;
				}
				k++;
			}

			if (!iterSwi)
				assert("LinearConjugateGradient function error : Exceed iteration number!");
		}

		static void LinearPreconditionedConjugateGradient(const MatrixX& A, const VectorX& b, const VectorX& x0, VectorX& x)
		{
			int n = A.rows();
			int maxIter = n * 2;
			int k = 0;
			bool iterSwi = false;

			MatrixX L(n, n), iL(n, n), iLt(n, n);
			VectorX r(n), y(n), y1(n), r1(n), p(n), Ap(n);
			Real alpha = 0, beta = 0;
			Real r1inner = 0;

			// initialization
			r = A * x0 - b;
			IncompleteCholeskyDecomposition(A, n, L);
			LowerMatrixInverse(L, n, iL);
			iLt = iL.transpose();
			y = iLt * iL * r; // 수정해야할듯
			p = -y; x = x0;

			while (k <= maxIter)
			{
				Ap = A * p;
				alpha = (r.transpose() * y);
				alpha /= (p.transpose()*Ap);
				x += alpha*p;
				r1 = r + alpha * Ap;
				y1 = iLt * iL * r1;
				r1inner = (r1.transpose()*y1);
				beta = r1inner / (r.transpose()*y);
				p = -y1 + beta*p;

				r = r1;
				y = y1;
				if ((r.transpose()*r) < 1E-10)
				{
					iterSwi = true;
					break;
				}
				k++;
			}

			if (!iterSwi)
				assert("LinearConjugateGradient function error : Exceed iteration number!");

		}
	}
}