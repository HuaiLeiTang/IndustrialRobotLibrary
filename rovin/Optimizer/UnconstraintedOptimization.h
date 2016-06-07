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

			// case 1
			//MatrixVectorMul(A, x0, n, n, r);
			//VectorSubtraction(r, b, n, r);
			
			// case 2
			Real r1inner = 0;
			r = A * x0 - b;

			// case 3
			//Real rinner = 0, r1inner = 0, pApinner = 0;
			//r.setZero();
			//for (int i = 0; i < n; i++)
			//{
			//	for (int j = 0; j < n; j++)
			//	{
			//		r(i) += A(i, j) * x0(j);
			//	}
			//	r(i) -= b(i);
			//}

			p = -r;
			x = x0;

			int maxIter = n * 2;

			while (k <= maxIter)
			{
				//std::cout << "k : " << k << std::endl;

				// case 1
				//MatrixVectorMul(A, p, n, n, Ap);
				//alpha = VectorInner(r, r, n) / VectorInner(p, Ap, n);
				//VectorAddition(x, p*alpha, n, x);
				//VectorAddition(r, Ap*alpha, n, r1);
				//beta = VectorInner(r1, r1, n) / VectorInner(r, r, n);
				//p *= beta;
				//VectorAddition(p, -r1, n, p);

				// case 2
				Ap = A * p;
				alpha = (r.transpose() * r);
				alpha /= (p.transpose()*Ap);
				x += alpha*p;
				r1 = r + alpha * Ap;
				r1inner = (r1.transpose()*r1);
				beta = r1inner / (r.transpose()*r);
				p = -r1 + beta*p;

				// case 3
				//rinner = 0;
				//r1inner = 0;
				//pApinner = 0;
				//Ap.setZero();

				//for (int i = 0; i < n; i++)
				//{
				//	for (int j = 0; j < n; j++)
				//	{
				//		Ap(i) += A(i, j) * p(j);
				//	}
				//}

				//for (int i = 0; i < n; i++)
				//	rinner += r(i) * r(i);

				//for (int i = 0; i < n; i++)
				//	pApinner += p(i) * Ap(i);

				//alpha = rinner / pApinner;

				//for (int i = 0; i < n; i++)
				//	x(i) = x(i) + alpha * p(i);

				//for (int i = 0; i < n; i++)
				//	r1(i) = r(i) + alpha * Ap(i);

				//for (int i = 0; i < n; i++)
				//	r1inner += r1(i) * r1(i);

				//beta = r1inner / rinner;
				//for (int i = 0; i < n; i++)
				//	p(i) = -r1(i) + beta * p(i);


				r = r1;
				//if (VectorInner(r, r, n) < 1E-15)
				//if((r.transpose()*r) < 1E-10)
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
	}
}