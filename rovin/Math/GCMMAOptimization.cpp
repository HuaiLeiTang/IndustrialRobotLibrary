#include "GCMMAOptimization.h"

#include "Common.h"
#include <iostream>
#include <time.h>

using namespace std;

namespace rovin
{
	void GCMMAOptimization::initialize(int xN, int ineqN)
	{
		setXN(xN);
		setIneqN(ineqN);

		setParameters(0.1, 0.5, 0.5, 0.7, 1.2);
		//setParameters(0.1, 0.5, 0.5, 0.8, 1.2);
		setCoefficients(1.0, VectorX(_ineqN).setZero(), VectorX(_ineqN).setConstant(20000.0), VectorX(_ineqN).setOnes());
		//setCoefficients(1.0, VectorX(_ineqN).setZero(), VectorX(_ineqN).setConstant(80000.0), VectorX(_ineqN).setOnes());

		_tolX = 1E-4;
		_tolFunc = 1E-4;
		_maxIterOL = 1000;
		_maxIterIL = 1000;
	}

	void GCMMAOptimization::setParameters(const Real & albefa, const Real & move0, const Real & asyinit, const Real & asydecr, const Real & asyincr)
	{
		_ALBEFA = albefa;
		_MOVE0 = move0;
		_ASYINIT = asyinit;
		_ASYDECR = asydecr;
		_ASYINCR = asyincr;
	}

	void GCMMAOptimization::setCoefficients(const Real & a0, const VectorX & ai, const VectorX & ci, const VectorX & di)
	{
		_a0 = a0;
		_ai = ai;
		_ci = ci;
		_di = di;
	}

	void GCMMAOptimization::setMinMax(const VectorX & minX, const VectorX & maxX)
	{
		_minX = minX;
		_maxX = maxX;
	}

	void GCMMAOptimization::solve(const VectorX & initialX)
	{

		// initialize
		// variables for outer loop

		// m1: values of 1 loop before, m2: values of 2 loops before
		VectorX xk(_xN), xkm1(_xN), xkm2(_xN);
		xk = initialX; xkm1 = xk; xkm2 = xk;

		//cout << "initialX" << endl << initialX << endl << endl;

		//VectorX low(_xN), lowm1(_xN), upp(_xN), uppm1(_xN), alpha(_xN), beta(_xN);
		_ollow.resize(_xN);
		_ollowm1.resize(_xN);
		_olupp.resize(_xN);
		_oluppm1.resize(_xN);
		_olalpha.resize(_xN);
		_olbeta.resize(_xN);

		// addition part
		tmpsubx.resize(_xN);
		tmpsublam.resize(_ineqN);
		tmpplam.resize(_xN);
		tmpqlam.resize(_xN);
		tmpdpsidx.resize(_xN);
		tmpgival.resize(_ineqN);

		G.resize(_ineqN, _xN);
		xadi.resize(_xN);
		bxdi.resize(_xN);
		ydi.resize(_ineqN);
		ldi.resize(_ineqN);

		Dx.resize(_xN, _xN); Dx.setZero();
		iDx.resize(_xN, _xN); iDx.setZero();
		Dly.resize(_ineqN, _ineqN); Dly.setZero();
		iDly.resize(_ineqN, _ineqN); iDly.setZero();
		iDy.resize(_ineqN);

		dxt.resize(_xN);
		dyt.resize(_ineqN);
		dlt.resize(_ineqN);
		dlyt.resize(_ineqN);

		DlyGDxGt.resize(_ineqN, _ineqN);
		iDlyGDxGt.resize(_ineqN, _ineqN);
		Lower.resize(_ineqN, _ineqN);
		iLower.resize(_ineqN, _ineqN);

		DxGtiDlyG.resize(_xN, _xN);
		iDxGtiDlyG.resize(_xN, _xN);
		Lower2.resize(_xN, _xN);
		iLower2.resize(_xN, _xN);

		MatSizeineqNbyxN.resize(_ineqN, _xN);
		MatSizexNbyineqN.resize(_xN, _ineqN);

		delx.resize(_xN);
		dellam.resize(_ineqN);

		// function evaluation variables
		VectorX f0val(1), f0valm1(1);		MatrixX df0dx(1, _xN), df0dxp(1, _xN), df0dxm(1, _xN); // originally, Real and VectorX respectively.
		VectorX fival(_ineqN);				MatrixX dfidx(_ineqN, _xN), dfidxp(_ineqN, _xN), dfidxm(_ineqN, _xN);
		f0valm1(0) = std::numeric_limits<Real>::max();

		// variables for inner loop
		VectorX xknu(_xN);
		// coefficients variables
		VectorX p0(_xN), q0(_xN);
		MatrixX pi(_ineqN, _xN), qi(_ineqN, _xN);
		Real r0, rho0;
		VectorX ri(_ineqN), rhoi(_ineqN);
		// function evaluation variables
		VectorX f0valknu(1), f0tvalknu(1), fivalknu(_ineqN), fitvalknu(_ineqN);

		int iterOL = 0, iterIL; // iter for outer/inner loop
		while (iterOL < _maxIterOL) // outer loop
		{
			//cout << "=== outer iter num: " << iterOL << endl;


			//cout << xk << endl << endl;
			//cout << xkm1 << endl << endl;
			//cout << xkm2 << endl << endl;

			calcLowUpp(iterOL, xk, xkm1, xkm2);
			//cout << "_olupp" << endl << _olupp << endl << endl;
			//cout << "_ollow" << endl << _ollow << endl << endl;
			calcAlphaBeta(xk);
			f0val = _objectFunc->func(xk);
			df0dx = _objectFunc->Jacobian(xk);
			fival = _ineqConstraint->func(xk);
			dfidx = _ineqConstraint->Jacobian(xk);
			calcPlusMinusMatrix(df0dx, df0dxp, df0dxm);
			calcPlusMinusMatrix(dfidx, dfidxp, dfidxm);


			//std::cout << low << std::endl << std::endl;
			//std::cout << upp << std::endl << std::endl;
			//std::cout << alpha << std::endl << std::endl;
			//std::cout << beta << std::endl << std::endl;

			//cout << f0val << endl;
			//cout << fival << endl;
			//cout << df0dx << endl << endl;
			//cout << df0dxp << endl << endl;
			//cout << df0dxm << endl << endl;
			//cout << dfidx << endl << endl;
			//cout << dfidxp << endl << endl;
			//cout << dfidxm << endl << endl;


			calcInitialRho(df0dx, dfidx, rho0, rhoi);

			//cout << rho0 << endl;
			//cout << rhoi << endl;

			//cout << "xk" << endl << xk << endl << endl;

			iterIL = 0;
			while (iterIL < _maxIterIL) // inner loop
			{
				//cout << "====== inner iter num: " << iterIL << endl;

				calcPQR(rho0, rhoi, df0dxp, df0dxm, dfidxp, dfidxm, xk, f0val, fival, p0, pi, q0, qi, r0, ri);

				//cout << p0 << endl << endl;
				//cout << pi << endl << endl;
				//cout << q0 << endl << endl;
				//cout << qi << endl << endl;
				//cout << r0 << endl << endl;
				//cout << ri << endl << endl;

				solveSubProblem(p0, pi, q0, qi, ri, xknu);

				//cout << xknu << endl << endl;


				if (testILSuccess(p0, pi, q0, qi, r0, ri, xknu, f0valknu, fivalknu, f0tvalknu, fitvalknu))
					break; // xknu is the optimal solution

							//cout << f0valknu << endl << endl;
							//cout << fivalknu << endl << endl;
							//cout << f0tvalknu << endl << endl;
							//cout << fitvalknu << endl << endl;

							//cout << rho0 << endl << endl;
							//cout << rhoi << endl << endl;

							// update rho0/rhoi
				updateRho0i(xknu, xk, f0valknu, fivalknu, f0tvalknu, fitvalknu, rho0, rhoi);

				//cout << rho0 << endl << endl;
				//cout << rhoi << endl << endl;

				iterIL++;
			}

			// terminate condition
			//cout << "abs(f0valm1(0) - f0valknu(0)) : " << abs(f0valm1(0) - f0valknu(0)) << endl;
			//cout << "(xkm1 - xknu).norm() : " << (xkm1 - xknu).norm() << endl;
			if (abs(f0valm1(0) - f0valknu(0)) < _tolFunc)
				break;
			if ((xkm1 - xknu).norm() < _tolX)
				break;

			// update
			xkm2 = xkm1;
			xkm1 = xk;
			xk = xknu;
			_ollowm1 = _ollow;
			_oluppm1 = _olupp;
			f0valm1(0) = f0valknu(0);

			//cout << "---------------" << endl;
			//cout << iterOL << endl << endl;
			//cout << xk << endl << endl;
			//cout << f0valknu << endl << endl;
			//cout << fivalknu << endl << endl << endl;


			iterOL++;
		}

		resultX = xknu;
		resultFunc = f0valknu(0);
	}


	void GCMMAOptimization::solveSubProblem(const VectorX & p0, const MatrixX & pi, const VectorX & q0, const MatrixX & qi, const VectorX & ri, VectorX & xout)
	{
		// input variable size:
		// VectorX p0(_xN), q0(_xN)
		// MatrixX pi(_ineqN, _xN), qi(_ineqN, _xN)
		// VectorX ri(_ineqN)
		// VectorX alpha(_xN), beta, low, upp

		// output variable size:
		// VectorX xout(_xN)


		// initialize
		VectorX bi = -ri;
		initializeSubProb();
		VectorX delw(_subDimW);
		Real tau;

		int iterSub = 0, maxIterSub = 1000;
		bool solFound = false;
		while (iterSub < maxIterSub)
		{
			//cout << "solveSubProblem iter num : " << iterSub << endl;
			//cout << "calcNormResidual_tmp(p0, pi, q0, qi, bi, delw, 0.0, -1) : " << calcNormResidual_tmp(p0, pi, q0, qi, bi, delw, 0.0, -1) << endl;
			//cout << "subx" << endl << _subx << endl << endl;

			// step 1: calculate gradient of w
			//calcGradientW(p0, pi, q0, qi, bi, delw);
			calcGradientW_tmp(p0, pi, q0, qi, bi, delw);

			//cout << "[delw]" << endl << delw << endl << endl;

			// step 2: calculate step length 'tau'
			tau = calcStepLength(delw, p0, pi, q0, qi, bi);

			//cout << "tau : " << tau << endl << endl;

			// step 3: update w
			_subw += tau * delw;
			separateFromW();

			// step 4: update epsilon
			//if (calcNormResidual(p0, pi, q0, qi, bi, delw, 0.0, -1) < 0.9 * _subeps)
			if (calcNormResidual_tmp(p0, pi, q0, qi, bi, delw, 0.0, -1) < 0.9 * _subeps)
				_subeps *= 0.1;

			// step 5: check terminate condition
			//if (_subeps <= 1E-6)
			if (_subeps <= 1E-4)
			{
				solFound = true;
				break;
			}

			iterSub++;
		}

		if (!solFound)
			LOG("exceeded max iteration number - 'solveSubProblem'");

		//cout << iterSub << endl;
		//cout << _subeps << endl;
		xout = _subx;


	}


	void GCMMAOptimization::calcLowUpp(int iter, const VectorX & xk, const VectorX & xkm1, const VectorX & xkm2)
	{
		if (iter == 0 || iter == 1) // first or second loop only
		{
			_ollow = xk - _ASYINIT * (_maxX - _minX);
			_olupp = xk + _ASYINIT * (_maxX - _minX);
		}
		else
		{
			for (int j = 0; j < _xN; j++)
			{
				if ((xk(j) - xkm1(j)) * (xkm1(j) - xkm2(j)) < 0)
				{
					_ollow(j) = xk(j) - _ASYDECR * (xkm1(j) - _ollowm1(j));
					_olupp(j) = xk(j) + _ASYDECR * (_oluppm1(j) - xkm1(j));
				}
				else if ((xk(j) - xkm1(j)) * (xkm1(j) - xkm2(j)) > 0)
				{
					_ollow(j) = xk(j) - _ASYINCR * (xkm1(j) - _ollowm1(j));
					_olupp(j) = xk(j) + _ASYINCR * (_oluppm1(j) - xkm1(j));
				}
				else
				{
					_ollow(j) = xk(j) - 1.0 * (xkm1(j) - _ollowm1(j));
					_olupp(j) = xk(j) + 1.0 * (_oluppm1(j) - xkm1(j));
				}
			}
		}
	}

	void GCMMAOptimization::calcAlphaBeta(const VectorX & xk)
	{
		Vector3 tmpVec;
		for (int j = 0; j < _xN; j++)
		{
			// for alpha
			tmpVec(0) = _minX(j);
			tmpVec(1) = _ollow(j) + _ALBEFA*(xk(j) - _ollow(j));
			tmpVec(2) = xk(j) - _MOVE0*(_maxX(j) - _minX(j));
			_olalpha(j) = tmpVec.maxCoeff();

			// for beta
			tmpVec(0) = _maxX(j);
			tmpVec(1) = _olupp(j) - _ALBEFA*(_olupp(j) - xk(j));
			tmpVec(2) = xk(j) + _MOVE0*(_maxX(j) - _minX(j));
			_olbeta(j) = tmpVec.minCoeff();
		}
	}

	void GCMMAOptimization::calcPlusMinusMatrix(const MatrixX & mat, MatrixX & matp, MatrixX & matm)
	{
		for (int i = 0; i < mat.rows(); i++)
		{
			for (int j = 0; j < mat.cols(); j++)
			{
				if (mat(i, j) > 0)
				{
					matp(i, j) = mat(i, j);
					matm(i, j) = 0;
				}
				else if (mat(i, j) < 0)
				{
					matp(i, j) = 0;
					matm(i, j) = -mat(i, j);
				}
				else
				{
					matp(i, j) = 0;
					matm(i, j) = 0;
				}
			}
		}
	}

	void GCMMAOptimization::calcInitialRho(const MatrixX & df0dx, const MatrixX & dfidx, Real & rho0, VectorX & rhoi)
	{
		// input variable size:	df0dx(1, _xN), dfidx(_ineqN, _xN)
		// output variable size: Real rho0, VectorX rhoi(_ineqN);

		Real tmpSum = 0;
		for (int j = 0; j < _xN; j++)
			tmpSum += abs(df0dx(0, j)) * (_maxX(j) - _minX(j));
		rho0 = 0.1 * tmpSum / (Real)_xN;

		for (int i = 0; i < _ineqN; i++)
		{
			tmpSum = 0;
			for (int j = 0; j < _xN; j++)
				tmpSum += abs(dfidx(i, j)) * (_maxX(j) - _minX(j));
			rhoi(i) = 0.1 * tmpSum / (Real)_xN;
			if (rhoi(i) < 1E-6)
				rhoi(i) = 1E-6;
		}
	}

	void GCMMAOptimization::calcPQR(const Real & rho0, const VectorX & rhoi, const MatrixX & df0dxp, const MatrixX & df0dxm, const MatrixX & dfidxp, const MatrixX & dfidxm, const VectorX & xk, const VectorX & f0val, const VectorX & fival, VectorX & p0, MatrixX & pi, VectorX & q0, MatrixX & qi, Real & r0, VectorX & ri)
	{
		// input variable size:
		// Real rho0
		// VectorX rhoi(_ineqN)
		// MatrixX df0dxp(1, _xN), df0dxm
		// MatrixX dfidxp(_ineqN, _xN), dfidxm
		// VectorX low(_xN), upp, xk
		// VectorX f0val(1), fival(_ineqN)

		// output variable size:
		// VectorX p0(_xN), q0(_xN)
		// MatrixX pi(_ineqN, _xN), qi(_ineqN, _xN)
		// Real r0,  VectorX ri(_ineqN)

		Real tmpSum = 0;

		for (int j = 0; j < _xN; j++)
		{
			p0(j) = (_olupp(j) - xk(j)) * (_olupp(j) - xk(j)) * (1.001 * df0dxp(0, j) + 0.001 * df0dxm(0, j) + rho0 / (_maxX(j) - _minX(j)));
			q0(j) = (xk(j) - _ollow(j)) * (xk(j) - _ollow(j)) * (0.001 * df0dxp(0, j) + 1.001 * df0dxm(0, j) + rho0 / (_maxX(j) - _minX(j)));

			tmpSum += p0(j) / (_olupp(j) - xk(j)) + q0(j) / (xk(j) - _ollow(j));
		}
		r0 = f0val(0) - tmpSum;

		for (int i = 0; i < _ineqN; i++)
		{
			tmpSum = 0;
			for (int j = 0; j < _xN; j++)
			{
				pi(i, j) = (_olupp(j) - xk(j)) * (_olupp(j) - xk(j)) * (1.001 * dfidxp(i, j) + 0.001 * dfidxm(i, j) + rhoi(i) / (_maxX(j) - _minX(j)));
				qi(i, j) = (xk(j) - _ollow(j)) * (xk(j) - _ollow(j)) * (0.001 * dfidxp(i, j) + 1.001 * dfidxm(i, j) + rhoi(i) / (_maxX(j) - _minX(j)));

				tmpSum += pi(i, j) / (_olupp(j) - xk(j)) + qi(i, j) / (xk(j) - _ollow(j));
			}
			ri(i) = fival(i) - tmpSum;
		}

	}

	void GCMMAOptimization::initializeSubProb()
	{
		// x, y, z, lam, xsi, eta, mu, zet, s (in order)
		_subDimW = _xN + _ineqN + 1 + _ineqN + _xN + _xN + _ineqN + 1 + _ineqN;
		_subw.resize(_subDimW);

		_subx.resize(_xN);
		_suby.resize(_ineqN);
		_sublam.resize(_ineqN);
		_subxsi.resize(_xN);
		_subeta.resize(_xN);
		_submu.resize(_ineqN);
		_subs.resize(_ineqN);

		resvec.resize(_subDimW);

		_subx = 0.5 * (_olalpha + _olbeta);
		_suby.setOnes();
		_subz = 1;
		_subzet = 1;
		_sublam.setOnes();
		_subs.setOnes();

		for (int j = 0; j < _xN; j++)
		{
			if (1 / (_subx(j) - _olalpha(j))>1)
				_subxsi(j) = 1 / (_subx(j) - _olalpha(j));
			else
			{
				_subxsi(j) = 1;
			}

			if (1 / (_olbeta(j) - _subx(j)) > 1)
				_subeta(j) = 1 / (_olbeta(j) - _subx(j));
			else
			{
				_subeta(j) = 1;
			}
		}

		for (int i = 0; i < _ineqN; i++)
		{
			if (0.5 * _ci(i) > 1)
				_submu(i) = 0.5 * _ci(i);
			else
				_submu(i) = 1;
		}

		combineToW(_subx, _suby, _subz, _sublam, _subxsi, _subeta, _submu, _subzet, _subs, _subw);



		_subeps = 1.0;
	}

	void GCMMAOptimization::separateFromW(void)
	{
		// save member variables x,y,z,lam,xsi... from member variable w
		int idx = 0;
		_subx = _subw.block(idx, 0, _xN, 1);			idx += _xN;
		_suby = _subw.block(idx, 0, _ineqN, 1);			idx += _ineqN;
		_subz = _subw(idx);								idx++;
		_sublam = _subw.block(idx, 0, _ineqN, 1);		idx += _ineqN;
		_subxsi = _subw.block(idx, 0, _xN, 1);			idx += _xN;
		_subeta = _subw.block(idx, 0, _xN, 1);			idx += _xN;
		_submu = _subw.block(idx, 0, _ineqN, 1);		idx += _ineqN;
		_subzet = _subw(idx);							idx++;
		_subs = _subw.block(idx, 0, _ineqN, 1);
	}

	void GCMMAOptimization::separateFromW(const VectorX & w, VectorX & x, VectorX & y, Real & z, VectorX & lam, VectorX & xsi, VectorX & eta, VectorX & mu, Real & zet, VectorX & s)
	{
		// save x,y,z,lam,xsi... from w
		int idx = 0;
		x = w.block(idx, 0, _xN, 1);			idx += _xN;
		y = w.block(idx, 0, _ineqN, 1);			idx += _ineqN;
		z = w(idx);								idx++;
		lam = w.block(idx, 0, _ineqN, 1);		idx += _ineqN;
		xsi = w.block(idx, 0, _xN, 1);			idx += _xN;
		eta = w.block(idx, 0, _xN, 1);			idx += _xN;
		mu = w.block(idx, 0, _ineqN, 1);		idx += _ineqN;
		zet = w(idx);							idx++;
		s = w.block(idx, 0, _ineqN, 1);
	}



	void GCMMAOptimization::combineToW(const VectorX & x, const VectorX & y, const Real & z, const VectorX & lam, const VectorX & xsi, const VectorX & eta, const VectorX & mu, const Real & zet, const VectorX & s, VectorX & w)
	{
		// save w from x,y,z,lam,xsi,...
		int idx = 0;
		w.block(idx, 0, _xN, 1) = x;			idx += _xN;
		w.block(idx, 0, _ineqN, 1) = y;			idx += _ineqN;
		w(idx) = z;								idx++;
		w.block(idx, 0, _ineqN, 1) = lam;		idx += _ineqN;
		w.block(idx, 0, _xN, 1) = xsi;			idx += _xN;
		w.block(idx, 0, _xN, 1) = eta;			idx += _xN;
		w.block(idx, 0, _ineqN, 1) = mu;		idx += _ineqN;
		w(idx) = zet;							idx++;
		w.block(idx, 0, _ineqN, 1) = s;
	}

	void GCMMAOptimization::calcGradientW(const VectorX & p0, const MatrixX & pi, const VectorX & q0, const MatrixX & qi, const VectorX & bi, VectorX & delw)
	{
		// output delw consists of delx, dely, delz, dellam, delxsi,...
		VectorX delx(_xN), dely(_ineqN), dellam(_ineqN), delxsi(_xN), deleta(_xN), delmu(_ineqN), dels(_ineqN);
		Real delz, delzet;

		MatrixX Psi(_xN, _xN), G(_ineqN, _xN);
		VectorX tmpplam(_xN), tmpqlam(_xN), tmpdpsidx(_xN), tmpgival(_ineqN);
		calcpqlam(p0, pi, _sublam, tmpplam);
		calcpqlam(q0, qi, _sublam, tmpqlam);
		calcdpsidx(tmpplam, tmpqlam, _subx, tmpdpsidx);
		calcgi(pi, qi, _subx, tmpgival);

		Psi.setZero();
		for (int j = 0; j < _xN; j++)
			Psi(j, j) = 2 * tmpplam(j) / pow(_olupp(j) - _subx(j), 3) + 2 * tmpqlam(j) / pow(_subx(j) - _ollow(j), 3);
		for (int i = 0; i < _ineqN; i++)
			for (int j = 0; j < _xN; j++)
				G(i, j) = pi(i, j) / pow(_olupp(j) - _subx(j), 2) - qi(i, j) / pow(_subx(j) - _ollow(j), 2);

		MatrixX xadi = (_subx - _olalpha).asDiagonal(); // <x-alpha>
		for (int i = 0; i < xadi.cols(); i++)
			xadi(i, i) = 1 / xadi(i, i);

		MatrixX bxdi = (_olbeta - _subx).asDiagonal().inverse(); // <beta -x>
		MatrixX ydi = _suby.asDiagonal().inverse(); // <y>
		MatrixX dd = _di.asDiagonal(); // <d>
		MatrixX ldi = _sublam.asDiagonal().inverse(); //<lamda>
		VectorX onesn(_xN);		onesn.setOnes();
		VectorX onesm(_ineqN);	onesm.setOnes();

		MatrixX Dx = Psi + xadi * _subxsi.asDiagonal() + bxdi * _subeta.asDiagonal();
		MatrixX Dy = dd + ydi * _submu.asDiagonal();
		MatrixX Dl = ldi * _subs.asDiagonal();
		MatrixX Dly = Dl + Dy.inverse();

		MatrixX iDx = MatrixX::Zero(Dx.rows(), Dx.cols());
		MatrixX iDy = MatrixX::Zero(Dy.rows(), Dy.cols());
		MatrixX iDly = MatrixX::Zero(Dly.rows(), Dly.cols());
		for (int i = 0; i < Dx.cols(); i++)
			iDx(i, i) = 1 / Dx(i, i);

		for (int i = 0; i < Dy.cols(); i++)
			iDy(i, i) = 1 / Dy(i, i);

		for (int i = 0; i < Dly.cols(); i++)
			iDly(i, i) = 1 / Dly(i, i);

		VectorX dxt = tmpdpsidx - _subeps * xadi * onesn + _subeps * bxdi * onesn;
		VectorX dyt = _ci + dd * _suby - _sublam - _subeps * ydi * onesm;
		Real dzt = _a0 - _sublam.transpose() * _ai - _subeps / _subz;
		VectorX dlt = tmpgival - _subz * _ai - _suby - bi + _subeps * ldi * onesm;
		VectorX dlyt = dlt + iDy * dyt;

		MatrixX smallA;
		VectorX smallx, smallb;
		if (_xN > _ineqN)
		{
			// equation (5.20)
			smallA.resize(_ineqN + 1, _ineqN + 1);
			smallx.resize(_ineqN + 1);
			smallb.resize(_ineqN + 1);
			smallA.block(0, 0, _ineqN, _ineqN) = Dly + G*iDx*G.transpose();
			smallA.block(0, _ineqN, _ineqN, 1) = _ai;
			smallA.block(_ineqN, 0, 1, _ineqN) = _ai.transpose();
			smallA(_ineqN, _ineqN) = -_subzet / _subz;
			smallb.block(0, 0, _ineqN, 1) = dlyt - G*iDx*dxt;
			smallb(_ineqN) = dzt;
			smallx = smallA.inverse() * smallb;
			dellam = smallx.block(0, 0, _ineqN, 1);
			delz = smallx(_ineqN);
			delx = -iDx*(G.transpose()*dellam + dxt);
		}
		else
		{
			// equation (5.22)
			smallA.resize(_xN + 1, _xN + 1);
			smallx.resize(_xN + 1);
			smallb.resize(_xN + 1);
			smallA.block(0, 0, _xN, _xN) = Dx + G.transpose()*iDly*G;
			smallA.block(0, _xN, _xN, 1) = -G.transpose()*iDly*_ai;
			smallA.block(_xN, 0, 1, _xN) = -_ai.transpose()*iDly*G;
			smallA(_xN, _xN) = _subzet / _subz + _ai.transpose()*iDly*_ai;
			smallb.block(0, 0, _xN, 1) = -dxt - G.transpose()*iDly*dlyt;
			smallb(_xN) = -dzt + _ai.transpose()*iDly*dlyt;
			smallx = smallA.inverse() * smallb;
			delx = smallx.block(0, 0, _xN, 1);
			delz = smallx(_xN);
			dellam = iDly*(G*delx - _ai*delz + dlyt);
		}
		dely = iDy * (dellam - dyt);
		delxsi = -xadi * (_subxsi.asDiagonal() * delx - _subeps * onesn) - _subxsi;
		deleta = bxdi * (_subeta.asDiagonal() * delx + _subeps * onesn) - _subeta;
		delmu = -ydi * (_submu.asDiagonal() * dely - _subeps * onesm) - _submu;
		delzet = -_subzet * delz / _subz - _subzet + _subeps / _subz;
		dels = -ldi * (_subs.asDiagonal() * dellam - _subeps * onesm) - _subs;

		combineToW(delx, dely, delz, dellam, delxsi, deleta, delmu, delzet, dels, delw);
	}

	void GCMMAOptimization::calcGradientW_tmp(const VectorX & p0, const MatrixX & pi, const VectorX & q0, const MatrixX & qi, const VectorX & bi, VectorX & delw)
	{
		// output delw consists of delx, dely, delz, dellam, delxsi,...
		delw.setZero();
		delx.setZero();
		dellam.setZero();

		calcpqlam(p0, pi, _sublam, tmpplam);
		calcpqlam(q0, qi, _sublam, tmpqlam);
		calcdpsidx(tmpplam, tmpqlam, _subx, tmpdpsidx);
		calcgi(pi, qi, _subx, tmpgival);

		//cout << "tmpdpsidx : " << tmpdpsidx << endl << endl;
		//cout << "tmpgival : " << tmpgival << endl << endl;
		//cout << "bi : " << bi << endl << endl;

		for (int i = 0; i < _ineqN; i++)
			for (int j = 0; j < _xN; j++)
				G(i, j) = pi(i, j) / pow(_olupp(j) - _subx(j), 2) - qi(i, j) / pow(_subx(j) - _ollow(j), 2);

		Gt = G.transpose();

		xadi = _subx - _olalpha;
		for (int i = 0; i < _xN; i++)
			xadi(i) = 1.0 / xadi(i);
		bxdi = _olbeta - _subx;
		for (int i = 0; i < _xN; i++)
			bxdi(i) = 1.0 / bxdi(i);
		for (int i = 0; i < _ineqN; i++)
			ydi(i) = 1.0 / _suby(i);
		for (int i = 0; i < _ineqN; i++)
			ldi(i) = 1.0 / _sublam(i);

		for (int i = 0; i < _xN; i++)
			Dx(i, i) = 2 * tmpplam(i) / pow(_olupp(i) - _subx(i), 3) + 2 * tmpqlam(i) / pow(_subx(i) - _ollow(i), 3) + xadi(i) * _subxsi(i) + bxdi(i) * _subeta(i);
		for (int i = 0; i < _xN; i++)
			iDx(i, i) = 1 / Dx(i, i);
		for (int i = 0; i < _ineqN; i++)
			iDy(i) = 1 / (_di(i) + ydi(i) * _submu(i));
		for (int i = 0; i < _ineqN; i++)
			Dly(i, i) = ldi(i) * _subs(i) + iDy(i);
		for (int i = 0; i < _ineqN; i++)
			iDly(i, i) = 1 / Dly(i, i);

		dxt = tmpdpsidx - _subeps * xadi + _subeps * bxdi;
		for (int i = 0; i < _ineqN; i++)
			dyt(i) = _ci(i) + _di(i) * _suby(i) - _sublam(i) - _subeps * ydi(i);
		dzt = _a0 - _sublam.transpose() * _ai - _subeps / _subz;
		dlt = tmpgival - _subz * _ai - _suby - bi + _subeps * ldi;
		for (int i = 0; i < _ineqN; i++)
			dlyt(i) = dlt(i) + iDy(i) * dyt(i);

		if (_xN > _ineqN)
		{
			// equation (5.20)
			iLower.setZero();
			iDlyGDxGt.setZero();

			// Calculate (Dly + G * iDx * G.transpose())
			Real sum = 0;
			for (int i = 0; i < _xN; i++)
			{
				for (int j = 0; j < _ineqN; j++)
					MatSizeineqNbyxN(j, i) = G(j, i) * iDx(i, i);
			}

			for (int i = 0; i < _ineqN; i++)
			{
				for (int j = 0; j < _ineqN; j++)
				{
					for (int k = 0; k < _xN; k++)
						sum += MatSizeineqNbyxN(i, k) * Gt(k, j);
					DlyGDxGt(i, j) = sum;
					sum = 0;
				}
				DlyGDxGt(i, i) += Dly(i, i);
			}

			// Calculate inverse of (Dly + G * iDx * G.transpose())
			for (int k = 0; k < _ineqN; k++) // Cholesky decomposition
			{
				Lower(k, k) = DlyGDxGt(k, k);
				for (int j = 0; j < k; j++)
					Lower(k, k) -= Lower(k, j) * Lower(k, j);
				Lower(k, k) = sqrt(Lower(k, k));

				for (int i = k + 1; i < _ineqN; i++)
				{
					Lower(i, k) = DlyGDxGt(i, k);
					for (int j = 0; j < k; j++)
						Lower(i, k) -= Lower(i, j) * Lower(k, j);
					Lower(i, k) /= Lower(k, k);
				}
			}

			for (int i = 0; i < _ineqN; i++) // Inverse of lower matrix
			{
				iLower(i, i) = 1 / Lower(i, i);
				for (int j = i + 1; j < _ineqN; j++)
				{
					for (int k = i; k < j; k++)
					{
						iLower(j, i) += Lower(j, k) * iLower(k, i);
					}
					iLower(j, i) /= -Lower(j, j);
				}
			}
			iLowert = iLower.transpose();

			for (int i = 0; i < _ineqN; i++) // multiple lower matrax & upper matrix
			{
				for (int j = i; j < _ineqN; j++)
				{
					for (int k = j; k < _ineqN; k++)
					{
						iDlyGDxGt(i, j) += iLowert(i, k)*iLower(k, j);
					}
					iDlyGDxGt(j, i) = iDlyGDxGt(i, j);
				}
			}

			// Calculate dellam, delz, delx
			for (int i = 0; i < _xN; i++) // dellam
			{
				for (int j = 0; j < _ineqN; j++)
					dellam(j) += G(j, i) * iDx(i, i) * dxt(i);
			}
			for (int i = 0; i < _ineqN; i++)
				dellam(i) = dlyt(i) - dellam(i);
			for (int i = 0; i < _ineqN; i++)
			{
				for (int j = 0; j < _ineqN; j++)
					delw(_xN + _ineqN + 1 + i) += iDlyGDxGt(i, j) * dellam(j);
			}
			delw(_xN + _ineqN) = -_subz * dzt / _subzet; // delz
			for (int i = 0; i < _xN; i++) // delx
			{
				for (int j = 0; j < _ineqN; j++)
					delx(i) += Gt(i, j) * delw(_xN + _ineqN + 1 + j);
			}
			for (int i = 0; i < _xN; i++)
				delx(i) += dxt(i);
			for (int i = 0; i < _xN; i++)
				delw(i) = -iDx(i, i) * delx(i);
		}
		else
		{
			// equation (5.22)
			iLower2.setZero();
			iDxGtiDlyG.setZero();

			// Calculate (Dx + G.transpose()*iDly*G) 
			Real sum = 0;
			for (int i = 0; i < _ineqN; i++)
			{
				for (int j = 0; j < _xN; j++)
					MatSizexNbyineqN(j, i) = Gt(j, i) * iDly(i, i);
			}
			for (int i = 0; i < _xN; i++)
			{
				for (int j = 0; j < _xN; j++)
				{
					for (int k = 0; k < _ineqN; k++)
						sum += MatSizexNbyineqN(i, k) * G(k, j);
					DxGtiDlyG(i, j) = sum;
					sum = 0;
				}
				DxGtiDlyG(i, i) += Dx(i, i);
			}

			// Calculate inverse of (Dx + G.transpose()*iDly*G)
			for (int k = 0; k < _xN; k++) // Cholesky decomposition
			{
				Lower2(k, k) = DxGtiDlyG(k, k);
				for (int j = 0; j < k; j++)
					Lower2(k, k) -= Lower2(k, j) * Lower2(k, j);
				Lower2(k, k) = sqrt(Lower2(k, k));

				for (int i = k + 1; i < _xN; i++)
				{
					Lower2(i, k) = DxGtiDlyG(i, k);
					for (int j = 0; j < k; j++)
						Lower2(i, k) -= Lower2(i, j) * Lower2(k, j);
					Lower2(i, k) /= Lower2(k, k);
				}
			}

			for (int i = 0; i < _xN; i++) // Inverse of lower matrix
			{
				iLower2(i, i) = 1 / Lower2(i, i);
				for (int j = i + 1; j < _xN; j++)
				{
					for (int k = i; k < j; k++)
					{
						iLower2(j, i) += Lower2(j, k) * iLower2(k, i);
					}
					iLower2(j, i) /= -Lower2(j, j);
				}
			}
			iLowert2 = iLower2.transpose();

			for (int i = 0; i < _xN; i++) // multiple lower matrax & upper matrix
			{
				for (int j = i; j < _xN; j++)
				{
					for (int k = j; k < _xN; k++)
					{
						iDxGtiDlyG(i, j) += iLowert2(i, k)*iLower2(k, j);
					}
					iDxGtiDlyG(j, i) = iDxGtiDlyG(i, j);
				}
			}

			// Calculate delx, delz, dellam
			for (int i = 0; i < _ineqN; i++) // delx
			{
				for (int j = 0; j < _xN; j++)
					delx(j) += Gt(j, i) * iDly(i, i) * dlyt(i);
			}
			for (int i = 0; i < _xN; i++)
				delx(i) = -dxt(i) - delx(i);
			for (int i = 0; i < _xN; i++)
			{
				for (int j = 0; j < _xN; j++)
					delw(i) += iDxGtiDlyG(i, j) * delx(j);
			}
			delw(_xN + _ineqN) = -_subz * dzt / _subzet; // delz
			for (int i = 0; i < _ineqN; i++) // dellam
			{
				for (int j = 0; j < _xN; j++)
					dellam(i) += G(i, j) * delw(j);
			}
			for (int i = 0; i < _ineqN; i++)
				dellam(i) += dlyt(i);
			for (int i = 0; i < _ineqN; i++)
				delw(_xN + _ineqN + 1 + i) = iDly(i, i) * dellam(i);
		}
		for (int i = 0; i < _ineqN; i++) // dely
			delw(i + _xN) = iDy(i) * (delw(i + _xN + _ineqN + 1) - dyt(i));
		for (int i = 0; i < _xN; i++) // delxsi
			delw(i + _xN + _ineqN + 1 + _ineqN) = -xadi(i) * (_subxsi(i) * delw(i) - _subeps) - _subxsi(i);
		for (int i = 0; i < _xN; i++) // deleta
			delw(i + _xN + _ineqN + 1 + _ineqN + _xN) = bxdi(i) * (_subeta(i) * delw(i) + _subeps) - _subeta(i);
		for (int i = 0; i < _ineqN; i++) // delmu
			delw(i + _xN + _ineqN + 1 + _ineqN + _xN + _xN) = -ydi(i) * (_submu(i) * delw(i + _xN) - _subeps) - _submu(i);
		delw(_xN + _ineqN + 1 + _ineqN + _xN + _xN + _ineqN) = -_subzet * delw(_xN + _ineqN) / _subz - _subzet + _subeps / _subz; // delzet
		for (int i = 0; i < _ineqN; i++) // dels
			delw(i + _xN + _ineqN + 1 + _ineqN + _xN + _xN + _ineqN + 1) = -ldi(i) * (_subs(i) * delw(i + _xN + _ineqN + 1) - _subeps) - _subs(i);
	}

	void GCMMAOptimization::calcpqlam(const VectorX & pq0, const MatrixX & pqi, const VectorX & lam, VectorX & pqlam)
	{
		// input variable size:
		// VectorX pq0: p0(_xN) or q0(_xN)
		// MatrixX pqi: pi(_ineqN, _xN) or qi(_ineqN, _xN)
		// VectorX lam(_ineqN)

		// output variable size:
		// VectorX pqlam: plam(_xN) or qlam(_xN)

		//cout << "lam" << endl << lam << endl << endl;
		//cout << "pq0" << endl << pq0 << endl << endl;
		//cout << "pqi" << endl << pqi << endl << endl;

		for (int j = 0; j < _xN; j++)
		{
			pqlam(j) = pq0(j);
			for (int i = 0; i < _ineqN; i++)
			{
				pqlam(j) += lam(i) * pqi(i, j);
			}
		}
	}

	void GCMMAOptimization::calcdpsidx(const VectorX & plam, const VectorX & qlam, const VectorX & x, VectorX & dpsidx)
	{
		//cout << "plam" << endl << plam << endl << endl;
		//cout << "qlam" << endl << qlam << endl << endl;
		//cout << "olupp" << endl << _olupp << endl << endl;
		//cout << "_ollow" << endl << _ollow << endl << endl;

		// equation 5.8
		for (int j = 0; j < _xN; j++)
			dpsidx(j) = plam(j) / pow(_olupp(j) - x(j), 2) - qlam(j) / pow(x(j) - _ollow(j), 2);
	}

	void GCMMAOptimization::calcgi(const MatrixX & pi, const MatrixX & qi, const VectorX & x, VectorX & gival)
	{
		// equation 5.2
		for (int i = 0; i < _ineqN; i++)
		{
			gival(i) = 0;
			for (int j = 0; j < _xN; j++)
				gival(i) += pi(i, j) / (_olupp(j) - x(j)) + qi(i, j) / (x(j) - _ollow(j));
		}
	}

	Real GCMMAOptimization::calcStepLength(const VectorX & delw, const VectorX & p0, const MatrixX & pi, const VectorX & q0, const MatrixX & qi, const VectorX & bi)
	{
		// minValLeq has the minimum value of sths (t <= sth)
		// maxValGeq has the maximum value of sths (t >= sth)
		Real minValLeq = 1, maxValGeq = -std::numeric_limits<Real>::max();
		Real coef = 0.99;

		// loop for 'xj + t*delxj - alphaj >= 0.01 * (xj - alphaj)
		for (int j = 0; j < _xN; j++)
		{
			if (delw(j) < 0)
			{
				// update minValLeq
				//if (coef*(_olalpha(j) - _subx(j)) / delw(j) < minValLeq)
				if ((std::abs(coef*(_olalpha(j) - _subx(j)) / delw(j)) > 0.00001) && coef*(_olalpha(j) - _subx(j)) / delw(j) < minValLeq)
					minValLeq = coef*(_olalpha(j) - _subx(j)) / delw(j);
			}
			else if (delw(j) > 0)
			{
				// update maxValGeq
				//if (coef*(_olalpha(j) - _subx(j)) / delw(j) > maxValGeq)
				if ((std::abs(coef*(_olalpha(j) - _subx(j)) / delw(j)) > 0.00001) && coef*(_olalpha(j) - _subx(j)) / delw(j) > maxValGeq)
					maxValGeq = coef*(_olalpha(j) - _subx(j)) / delw(j);
			}
		}

		// loop for 'betaj - (xj + t*delxj) >= 0.01 * (betaj - xj)'
		for (int j = 0; j < _xN; j++)
		{
			if (delw(j) > 0)
			{
				// update minValLeq
				//if (coef*(_olbeta(j) - _subx(j)) / delw(j) < minValLeq)
				if ((std::abs(coef*(_olbeta(j) - _subx(j)) / delw(j)) > 0.00001) && coef*(_olbeta(j) - _subx(j)) / delw(j) < minValLeq)
					minValLeq = coef*(_olbeta(j) - _subx(j)) / delw(j);
			}
			else if (delw(j) < 0)
			{
				// update maxValGeq
				//if (coef*(_olbeta(j) - _subx(j)) / delw(j) > maxValGeq)
				if ((std::abs(coef*(_olbeta(j) - _subx(j)) / delw(j)) > 0.00001) && coef*(_olbeta(j) - _subx(j)) / delw(j) > maxValGeq)
					maxValGeq = coef*(_olbeta(j) - _subx(j)) / delw(j);
			}
		}

		// loop for (y,z,lam,...) + t * (dely, delz, dellam, ...) >= 0.01 * (y, z, lam, ...)
		for (int idx = _xN; idx < _subDimW; idx++)
		{
			if (delw(idx) < 0)
			{
				// update minValLeq
				//if (-coef * _subw(idx) / delw(idx) < minValLeq)
				if ((std::abs(-coef * _subw(idx) / delw(idx)) > 0.00001) && -coef * _subw(idx) / delw(idx) < minValLeq)
					minValLeq = -coef * _subw(idx) / delw(idx);
			}
			else if (delw(idx) > 0)
			{
				// update maxValGeq
				//if (-coef * _subw(idx) / delw(idx) > maxValGeq)
				if ((std::abs(-coef * _subw(idx) / delw(idx)) > 0.00001) && -coef * _subw(idx) / delw(idx) > maxValGeq)
					maxValGeq = -coef * _subw(idx) / delw(idx);
			}
		}

		LOGIF(maxValGeq <= minValLeq, "'minValLeq' must be larger than 'maxValGeq'");

		//cout << "maxValGeq : " << maxValGeq << endl;
		//cout << "minVelLeg : " << minValLeq << endl;

		//Real resNorm = calcNormResidual(p0, pi, q0, qi, bi, delw, 0.0, 2);
		Real resNorm = calcNormResidual_tmp(p0, pi, q0, qi, bi, delw, 0.0, 2);
		bool tauFound = false;
		int iterTau = 0, maxIterTau = 1000;

		//cout << "minValLeq : " << minValLeq << endl;
		while (iterTau < maxIterTau)
		{
			//if (calcNormResidual(p0, pi, q0, qi, bi, delw, minValLeq, 2) < resNorm)
			if (calcNormResidual_tmp(p0, pi, q0, qi, bi, delw, minValLeq, 2) < resNorm)
			{
				//cout << "minValLeg : " << minValLeq << endl;
				//cout << "calcNormResidual_tmp(p0, pi, q0, qi, bi, delw, minValLeq, 2) : " << calcNormResidual_tmp(p0, pi, q0, qi, bi, delw, minValLeq, 2) << endl;
				//cout << "calcNormResidual_tmp(p0, pi, q0, qi, bi, delw, 0.0, 2) : " << calcNormResidual_tmp(p0, pi, q0, qi, bi, delw, 0.0, 2) << endl;

				tauFound = true;
				break;
			}
			minValLeq *= 0.5;
			iterTau++;

			if (iterTau == maxIterTau)
			{
				//cout << "init : " << init << endl;
				//cout << "maxValGeq : " << maxValGeq << endl;
				//cout << "minValLeg : " << minValLeq << endl;
				//cout << calcNormResidual_tmp(p0, pi, q0, qi, bi, delw, minValLeq, 2) << endl;
				//cout << calcNormResidual_tmp(p0, pi, q0, qi, bi, delw, 0.0, 2) << endl;
				//cout << endl;
			}

		}
		LOGIF(tauFound, "failed to find step length 'tau'");

		return minValLeq;
	}

	Real GCMMAOptimization::calcNormResidual(const VectorX & p0, const MatrixX & pi, const VectorX & q0, const MatrixX & qi, const VectorX & bi, const VectorX & delw, Real stepLength, int normCh)
	{
		// output residual vector consists of deltax, deltay, deltaz, deltalam, deltaxsi,...
		VectorX resvec(_subDimW);
		VectorX deltax(_xN), deltay(_ineqN), deltalam(_ineqN), deltaxsi(_xN), deltaeta(_xN), deltamu(_ineqN), deltas(_ineqN); Real deltaz, deltazet;

		// del: separate variable of delw
		VectorX delx(_xN), dely(_ineqN), dellam(_ineqN), delxsi(_xN), deleta(_xN), delmu(_ineqN), dels(_ineqN); Real delz, delzet;
		separateFromW(delw, delx, dely, delz, dellam, delxsi, deleta, delmu, delzet, dels);


		// tmpsub = _sub + stepLength * del
		VectorX tmpsubx(_xN), tmpsuby(_ineqN), tmpsublam(_ineqN), tmpsubxsi(_xN), tmpsubeta(_xN), tmpsubmu(_ineqN), tmpsubs(_ineqN); Real tmpsubz, tmpsubzet;
		tmpsubx = _subx + stepLength * delx;
		tmpsuby = _suby + stepLength * dely;
		tmpsubz = _subz + stepLength * delz;
		tmpsublam = _sublam + stepLength * dellam;
		tmpsubxsi = _subxsi + stepLength * delxsi;
		tmpsubeta = _subeta + stepLength * deleta;
		tmpsubmu = _submu + stepLength * delmu;
		tmpsubzet = _subzet + stepLength * delzet;
		tmpsubs = _subs + stepLength * dels;

		// dpsidx, gival 다 새로 구해야함, tmpsub variable 로.....
		VectorX tmpplam(_xN), tmpqlam(_xN), tmpdpsidx(_xN), tmpgival(_ineqN);
		calcpqlam(p0, pi, tmpsublam, tmpplam);
		calcpqlam(q0, qi, tmpsublam, tmpqlam);
		calcdpsidx(tmpplam, tmpqlam, tmpsubx, tmpdpsidx);
		calcgi(pi, qi, tmpsubx, tmpgival);


		deltax = tmpdpsidx - tmpsubxsi + tmpsubeta;
		for (int i = 0; i < _ineqN; i++)
		{
			deltay(i) = _ci(i) + _di(i) * tmpsuby(i) - tmpsublam(i) - tmpsubmu(i);
			deltamu(i) = tmpsubmu(i) * tmpsuby(i) - _subeps;
			deltas(i) = tmpsublam(i) * tmpsubs(i) - _subeps;
		}
		deltaz = _a0 - tmpsubzet - tmpsublam.transpose() * _ai;
		deltalam = tmpgival - tmpsubz * _ai - tmpsuby + tmpsubs - bi;
		for (int j = 0; j < _xN; j++)
		{
			deltaxsi(j) = tmpsubxsi(j) * (tmpsubx(j) - _olalpha(j)) - _subeps;
			deltaeta(j) = tmpsubeta(j) * (_olbeta(j) - tmpsubx(j)) - _subeps;
		}
		deltazet = tmpsubzet * tmpsubz - _subeps;


		combineToW(deltax, deltay, deltaz, deltalam, deltaxsi, deltaeta, deltamu, deltazet, deltas, resvec);

		Real normVal;
		switch (normCh)
		{
		case 2: // Euclidean norm
			normVal = resvec.norm();
			break;
		case -1: // infinite norm
			normVal = resvec.cwiseAbs().maxCoeff();
			break;
		default:
			LOG("wrong choice - function 'calcNormResidual'");
			break;
		}
		return normVal;
	}

	Real GCMMAOptimization::calcNormResidual_tmp(const VectorX& p0, const MatrixX& pi, const VectorX& q0, const MatrixX& qi, const VectorX& bi, const VectorX& delw, Real stepLength, int normCh)
	{
		//VectorX resvec(_subDimW);

		//cout << "subx" << endl << _subx << endl << endl;

		//VectorX tmpsubx(_xN), tmpsublam(_ineqN);
		for (int i = 0; i < _xN; i++)
			tmpsubx(i) = _subx(i) + stepLength * delw(i);
		for (int i = 0; i < _ineqN; i++)
			tmpsublam(i) = _sublam(i) + stepLength * delw(i + _xN + _ineqN + 1);

		//cout << "tmpsubx" << endl << tmpsubx << endl << endl;

		//VectorX tmpplam(_xN), tmpqlam(_xN), tmpdpsidx(_xN), tmpgival(_ineqN);
		calcpqlam(p0, pi, tmpsublam, tmpplam);
		calcpqlam(q0, qi, tmpsublam, tmpqlam);
		calcdpsidx(tmpplam, tmpqlam, tmpsubx, tmpdpsidx);
		calcgi(pi, qi, tmpsubx, tmpgival);

		//cout << "tmpdpsidx" << endl << tmpdpsidx << endl << endl;
		//cout << "_subxsi" << endl << _subxsi << endl;
		//cout << "_subeta" << endl << _subeta << endl;

		int idx = 0;
		for (int i = 0; i < _xN; i++)
			resvec(i + idx) = tmpdpsidx(i) - _subxsi(i) + stepLength * delw(i + _xN + _ineqN + 1 + _ineqN) + _subeta(i) + stepLength * delw(i + _xN + _ineqN + 1 + _ineqN + _xN);
		idx += _xN;
		for (int i = 0; i < _ineqN; i++)
			resvec(i + idx) = _ci(i) + _di(i) * (_suby(i) + stepLength * delw(i + _xN)) - tmpsublam(i) - (_submu(i) + stepLength * delw(i + _xN + _ineqN + 1 + _ineqN + _xN + _xN));
		idx += _ineqN;
		resvec(idx) = _a0 - (_subzet + stepLength * delw(_xN + _ineqN + 1 + _ineqN + _xN + _xN + _ineqN)) - tmpsublam.transpose() * _ai;
		idx += 1;
		for (int i = 0; i < _ineqN; i++)
			resvec(i + idx) = tmpgival(i) - (_subz + stepLength * delw(_xN + _ineqN)) * _ai(i) - (_suby(i) + stepLength * delw(i + _xN)) + (_subs(i) + stepLength * delw(i + _xN + _ineqN + 1 + _ineqN + _xN + _xN + _ineqN + 1)) - bi(i);
		idx += _ineqN;
		for (int i = 0; i < _xN; i++)
			resvec(i + idx) = (_subxsi(i) + stepLength * delw(i + _xN + _ineqN + 1 + _ineqN)) * (tmpsubx(i) - _olalpha(i)) - _subeps;
		idx += _xN;
		for (int i = 0; i < _xN; i++)
			resvec(i + idx) = (_subeta(i) + stepLength * delw(i + _xN + _ineqN + 1 + _ineqN + _xN)) * (_olbeta(i) - tmpsubx(i)) - _subeps;
		idx += _xN;
		for (int i = 0; i < _ineqN; i++)
			resvec(i + idx) = (_submu(i) + stepLength * delw(i + _xN + _ineqN + 1 + _ineqN + _xN + _xN)) * (_suby(i) + stepLength * delw(i + _xN)) - _subeps;
		idx += _ineqN;
		resvec(idx) = (_subzet + stepLength * delw(_xN + _ineqN + 1 + _ineqN + _xN + _xN + _ineqN)) * (_subz + stepLength * delw(_xN + _ineqN)) - _subeps;
		idx += 1;
		for (int i = 0; i < _ineqN; i++)
			resvec(i + idx) = tmpsublam(i) * (_subs(i) + stepLength * delw(i + _xN + _ineqN + 1 + _ineqN + _xN + _xN + _ineqN + 1)) - _subeps;

		//cout << "resvec" << endl << resvec << endl << endl;

		Real normVal;
		switch (normCh)
		{
		case 2: // Euclidean norm
			normVal = resvec.norm();
			break;
		case -1: // infinite norm
			normVal = resvec.cwiseAbs().maxCoeff();
			break;
		default:
			LOG("wrong choice - function 'calcNormResidual'");
			break;
		}
		return normVal;

		return 0;
	}

	bool GCMMAOptimization::testILSuccess(const VectorX & p0, const MatrixX & pi, const VectorX & q0, const MatrixX & qi, const Real & r0, const VectorX & ri, const VectorX & testx, VectorX & f0valknu, VectorX & fivalknu, VectorX & f0tvalknu, VectorX & fitvalknu)
	{
		//VectorX f0val(1), f0tval(1), fival(_ineqN), fitval(_ineqN);
		f0valknu = _objectFunc->func(testx);
		fivalknu = _ineqConstraint->func(testx);
		calcf0tilde(p0, q0, r0, testx, f0tvalknu);
		calcfitilde(pi, qi, ri, testx, fitvalknu);

		bool ret = false;

		if (f0tvalknu(0) >= f0valknu(0))
		{
			ret = true;
			for (int i = 0; i < _ineqN; i++)
			{
				if (fitvalknu(i) < fivalknu(i))
				{
					ret = false;
					break;
				}
			}
		}
		return ret;
	}


	void GCMMAOptimization::calcf0tilde(const VectorX & p0, const VectorX & q0, const Real & r0, const VectorX & x, VectorX & f0tval)
	{
		f0tval(0) = r0;
		for (int j = 0; j < _xN; j++)
			f0tval(0) += p0(j) / (_olupp(j) - x(j)) + q0(j) / (x(j) - _ollow(j));
	}

	void GCMMAOptimization::calcfitilde(const MatrixX & pi, const MatrixX & qi, const VectorX & ri, const VectorX & x, VectorX & fitval)
	{
		for (int i = 0; i < _ineqN; i++)
		{
			fitval(i) = ri(i);
			for (int j = 0; j < _xN; j++)
				fitval(i) += pi(i, j) / (_olupp(j) - x(j)) + qi(i, j) / (x(j) - _ollow(j));
		}
	}

	void GCMMAOptimization::updateRho0i(const VectorX & xknu, const VectorX & xk, const VectorX & f0valknu, const VectorX & fivalknu, const VectorX & f0tvalknu, const VectorX & fitvalknu, Real & rho0, VectorX & rhoi)
	{
		Real dknu = 0;
		for (int j = 0; j < _xN; j++)
			dknu += (_olupp(j) - _ollow(j)) * (xknu(j) - xk(j)) * (xknu(j) - xk(j)) / ((_olupp(j) - xknu(j))*(xknu(j) - _ollow(j))*(_maxX(j) - _minX(j)));

		Real deltaknu0 = (f0valknu(0) - f0tvalknu(0)) / dknu;
		if (deltaknu0 > 0)
			if (1.1*(rho0 + deltaknu0) < 10 * rho0)
				rho0 = 1.1*(rho0 + deltaknu0);
			else
				rho0 = 10 * rho0;
		VectorX deltaknui = (fivalknu - fitvalknu) / dknu;
		for (int i = 0; i < _ineqN; i++)
		{
			if (deltaknui(i) > 0)
				if (1.1*(rhoi(i) + deltaknui(i) < 10 * rhoi(i)))
					rhoi(i) = 1.1*(rhoi(i) + deltaknui(i));
				else
					rhoi(i) = 10 * rhoi(i);
		}
	}

}