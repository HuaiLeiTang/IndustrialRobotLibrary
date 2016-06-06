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
			int k = 0;

			VectorX r(n), r1(n), p(n), Ap(n);
			Real alpha = 0, beta = 0;
			bool iterSwi = false;

			MatrixVectorMul(A, x0, n, n, r);
			VectorSubtraction(r, b, n, r);
			p = -r;
			x = x0;

			while (k <= n)
			{
				MatrixVectorMul(A, p, n, n, Ap);
				alpha = VectorInner(r, r, n) / VectorInner(p, Ap, n);
				VectorAddition(x, p*alpha, n, x);
				VectorAddition(r, Ap*alpha, n, r1);
				beta = VectorInner(r1, r1, n) / VectorInner(r, r, n);
				p *= beta;
				VectorAddition(p, -r1, n, p);

				r = r1;


				if (VectorInner(r, r, n) < 1E-15)
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