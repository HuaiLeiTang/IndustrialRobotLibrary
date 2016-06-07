#include "convexQP.h"

#define EPS_COMPARE 1E-4

#include <iostream>
using namespace std;

namespace rovin
{
	namespace irOpt
	{



		bool convexQP::solveInequalityConstrainedQP(void)
		{
			//checkDimension();

			initialize();

			// outer loop
			int iterOL = 0, maxIterOL = 1000, iterIL, maxIterIL = 1000;
			while (iterOL++ < maxIterOL) 
			{
				// check solution found
				if (checkSuccess())
					break;

				// find index 'r': index of value of inequality has least negative value in idxWSdual
				selectIndexR();
				//cout << _idxr << endl << endl;

				// inner loop
				iterIL = 0;
				while ((iterIL++ < maxIterIL) && RealLess(_ineqCon(_idxr), 0.0, EPS_COMPARE))
				{
					// calculate p and v
					calcPandV();
					//cout << _p << endl << endl;
					//cout << _v << endl << endl;					

					if (RealEqual(VectorInner(_A.row(_idxr), _p, _xN), 0.0, EPS_COMPARE))
					{
						if (checkFailure())
						{
							return false; // failure
						}
						else
						{
							// compute step length _t1, remove constraint _idx1							
							calcTandIndex();
							updateMu(_t1);
							removeConFromWS(_idx1);

						}
					}
					else
					{
						// compute step length tD (_t1) and tP (_t2)
						calcTandIndex();
						_t1 = Min(_t1, RealMax);
						_t2 = -_ineqCon(_idxr) / VectorInner(_A.row(_idxr), _p, _xN);
						_t3 = Min(_t1, _t2);
						
						for (int i = 0; i < _xN; i++)
							_x(i) += _t3 * _p(i);
						updateMu(_t3);

						if (RealLessEqual(_t2, _t1, EPS_COMPARE))
							appendConToWS(_idxr);
						else
							removeConFromWS(_idx1);

					}


					// update
					_ineqCon = _A * _x - _b;
				}
			}




			_resultX = _x;
			_resultmu = _mu;


			return true;
		}

		void convexQP::initialize(void)
		{
			// make initial guess
			//_x.resize(_xN);
			//MatrixX iG(_xN, _xN);
			//PosDefMatrixInverse(_G, _xN, iG);
			//MatrixVectorMul(iG, _c, _xN, _xN, _x);
			//for (int i = 0; i < _xN; i++)
			//	_x(i) *= -1;
			_x = -_G.inverse() * _c;


			// initial mu (Lagrange multiplier)
			_mu.setZero(_ineqN);


			// initial working set
			_noWSdual = _ineqN;
			_noWS = 0;
			_idxWSdual.clear();
			_idxWS.clear();
			for (int i = 0; i < _ineqN; i++)
				_idxWSdual.push_back(i);

			// initial values of constraints
			_ineqCon = _A * _x - _b;


			// 다른 멤버변수들 여기서 초기화 하기!
			// G, c, A, b, _xN, _ineqN 말고
		}

		bool convexQP::checkSuccess(void)
		{
			for (int i = 0; i < _ineqN; i++) {
				if (RealLess(_ineqCon(i), 0.0, EPS_COMPARE)) {
					return false;
				}
			}
			return true;
		}

		bool convexQP::checkFailure(void)
		{
			for (int i = 0; i < _noWS; i++) {
				if (RealLess(_v(i), 0.0, EPS_COMPARE)) {
					return false;
				}
			}
			return true; // failure of convex QP algorithm
		}

		void convexQP::selectIndexR(void)
		{
			bool set = false;
			for (_iterList = _idxWSdual.begin(); _iterList != _idxWSdual.end(); _iterList++)
			{
				if (_ineqCon(*_iterList) < 0)
				{
					if (!set)
					{
						_idxr = *_iterList;
						set = true;
					}
					else
					{
						if (_ineqCon(*_iterList) > _ineqCon(_idxr))
							_idxr = *_iterList;
					}
				}
			}
		}

		void convexQP::calcPandV(void)
		{
			// 여기 좀 빠르게 고쳐야 함.. 일단 아이겐 라이브러리 inverse....
			int i, j;

			_p.resize(_xN);
			if (_noWS == 0)
			{
				_v.resize(0);
				_tmpA = _G;
				_tmpb = _A.row(_idxr);
				_p = _tmpA.inverse() * _tmpb;
			}
			else
			{
				_v.resize(_noWS);
				_tmpA.resize(_xN + _noWS, _xN + _noWS);
				_tmpA.setZero();
				_tmpb.resize(_xN + _noWS);
				_tmpb.setZero();

				for (i = 0; i < _xN; i++) {
					for (j = 0; j < _xN; j++) {
						_tmpA(i, j) = _G(i, j);
					}
				}
				i = 0;
				for (_iterList = _idxWS.begin(); _iterList != _idxWS.end(); _iterList++) {
					for (j = 0; j < _xN; j++) {
						_tmpA(_xN + i, j) = -_A(*_iterList, j);
						_tmpA(j, _xN + i) = -_A(*_iterList, j);
					}
					i++;
				}

				for (i = 0; i < _xN; i++) {
					_tmpb(i) = _A(_idxr, i);
				}
				// 여기 들어올때 A,b확인하기...
				//cout<<_tmpA<<endl<<endl;
				//cout<<_tmpb<<endl<<endl;

				_tmpx = _tmpA.inverse() * _tmpb;
				for (i = 0; i < _xN; i++)
					_p(i) = _tmpx(i);
				for (i = 0; i < _noWS; i++)
					_v(i) = _tmpx(_xN + i);
			}
		}

		void convexQP::calcTandIndex(void)
		{
			// calculate _t1, _idx1
			bool set = false;
			int i = 0;
			for (_iterList = _idxWS.begin(); _iterList != _idxWS.end(); _iterList++)
			{
				if (RealLess(_v(i), 0.0, EPS_COMPARE))
				{
					if (!set)
					{
						_idx1 = *_iterList;
						_t1 = -_mu(_idx1) / _v(i);
						set = true;
					}
					else
					{
						if (-_mu(*_iterList) / _v(i) < _t1)
						{
							_idx1 = *_iterList;
							_t1 = -_mu(_idx1) / _v(i);
						}
					}
				}
				i++;
			}

			if (!set)
			{
				_t1 = RealMax;
				_idx1 = -1;
			}
		}

		void convexQP::updateMu(const Real & t)
		{
			// update _mu with t
			int i = 0;
			for (_iterList = _idxWS.begin(); _iterList != _idxWS.end(); _iterList++)
			{
				_mu(*_iterList) += t * _v(i);
				i++;
			}
			_mu(_idxr) += t;
		}

		void convexQP::appendConToWS(int idx)
		{
			// append idx to _idxWS, remove idx from _idxWSdual

			for (_iterList = _idxWS.begin(); _iterList != _idxWS.end(); _iterList++)
			{
				if (*_iterList > idx)
					break;
			}
			_idxWS.insert(_iterList, idx);
			_noWS++;

			_iterList = find(_idxWSdual.begin(), _idxWSdual.end(), idx);
			_idxWSdual.erase(_iterList);
			_noWSdual--;

		}

		void convexQP::removeConFromWS(int idx)
		{
			// remove idx from _idxWS, append idx to _idxWSdual

			_iterList = find(_idxWS.begin(), _idxWS.end(), idx);
			_idxWS.erase(_iterList);
			_noWS--;

			for (_iterList = _idxWSdual.begin(); _iterList != _idxWSdual.end(); _iterList++)
			{
				if (*_iterList > idx)
					break;
			}
			_idxWSdual.insert(_iterList, idx);
			_noWSdual++;
		}



	}
}
