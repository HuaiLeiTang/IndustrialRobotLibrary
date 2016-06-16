#pragma once

#include "rovin\Math\Constant.h"
#include "rovin\Math\Common.h"
#include <iostream>

#include "rovin\Math\Function.h"

using namespace rovin;
using namespace std;

namespace irLib
{
	namespace Opt
	{
		class NewtonRaphson;
		class NewtonRaphsonfcn;

		typedef void(*fcn)(const VectorX& x, VectorX& f, MatrixX& ig, void* f_data);

		class NewtonRaphsonfcn
		{
		private:
			int _xN;
			int _fN;
			Real _tol;
			fcn _fcn;
			void* _f_data;
			VectorX resultX;
		public:
			NewtonRaphsonfcn(int xN, int fN, Real tol = 1E-10) : _xN(xN), _fN(fN), _tol(tol)
			{
				resultX.resize(xN);
				resultX.setZero();
			}
			void setxN(const int xN) { _xN = xN; }
			void setfN(const int fN) { _fN = fN; }
			void setfunction(fcn fcnt) { _fcn = fcnt; }
			void setfdata(void* f_data) { _f_data = f_data; }
			void settolerance(const Real tol) { _tol = tol; }

			void solve(const VectorX& initX, int maxIter = 1000)
			{
				int iter = 0;
				VectorX f(_fN);
				MatrixX g(_fN, _xN);
				
				VectorX tmp(_fN);

				resultX = initX;

				//cout << "xN : " << _xN << '\t' << "fN : " << _fN << endl;

				while (iter < maxIter)
				{
					_fcn(resultX, f, g, _f_data);
					tmp.setZero();
					for (int i = 0; i < _fN; i++)
					{
						for (int j = 0; j < _xN; j++)
						{
							//cout << "i : " << i << '\t' << "j : " << j << endl;
							tmp(i) += g(i, j) * f(j);
						}
					}
					resultX -= tmp;

					if (f.norm() < _tol)
						break;
					iter++;
				}
				cout << iter << endl;
			}
			const VectorX& getResultX() const { return resultX; }

		};

		class NewtonRaphson
		{
		private:
			int _xN;
			int _fN;
			Real _tol;
			FunctionPtr _fcn;
			VectorX resultX;
		public:
			NewtonRaphson(int xN, int fN, Real tol = 1E-10) : _xN(xN), _fN(fN), _tol(tol) 
			{
				resultX.resize(xN);
				resultX.setZero();
			}

			void setxN(const int xN) { _xN = xN; }
			void setfN(const int fN) { _fN = fN; }
			void setfunction(FunctionPtr fcn) { _fcn = fcn; }
			void settolerance(const Real tol) { _tol = tol; }

			void solve(const VectorX& initX, int maxIter = 1000)
			{
				int iter = 0;
				VectorX f(_fN);

				MatrixX g(_fN, _xN);
				VectorX tmp(_fN);

				resultX = initX;

				//cout << "xN : " << _xN << '\t' << "fN : " << _fN << endl;
				
				while (iter < maxIter)
				{	
					f = _fcn->func(resultX);
					
					//resultX -= _fcn->InverseJacobian(resultX) * f;
					g = _fcn->InverseJacobian(resultX);
					tmp.setZero();
					for (int i = 0; i < _fN; i++)
					{
						for (int j = 0; j < _xN; j++)
						{
							//cout << "i : " << i << '\t' << "j : " << j << endl;
							tmp(i) += g(i, j) * f(j);
						}
					}
					resultX -= tmp;

					if (f.norm() < _tol)
						break;
					iter++;
				}
				cout << iter << endl;
			}
			const VectorX& getResultX() const { return resultX; }
		};

		static void LinearConjugateGradient(const MatrixX& A, const VectorX& b, const VectorX& x0, VectorX& x)
		{
			int n = A.rows();
			int nn = 0;
			int maxIter = n * 10;
			int k = 0;
			bool iterSwi = false;

			VectorX r(n), r1(n), p0(n), p(n), Ap(n);
			Real alpha = 0, beta = 0;
			Real r1inner = 0;
			
			// initialization
			r = A * x0 - b;
			p0 = -r; p = -r; x = x0;

			while (k < maxIter)
			{
				cout << "k : " << k << endl;
				Ap = A * p;
				alpha = (r.transpose() * r);
				alpha /= (p.transpose()*Ap);
				x += alpha*p;
				r1 = r + alpha * Ap;
				r1inner = (r1.transpose()*r1);
				beta = r1inner / (r.transpose()*r);
				p = -r1 + beta*p;

				r = r1;
				cout << "r1inner : " << r1inner << endl;
				if(r1inner < 1E-10)
				{
					iterSwi = true;
					break;
				}

				//if (nn == (n - 1))
				//{
				//	p = p0;
				//	nn = 0;
				//}
				//nn++;

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