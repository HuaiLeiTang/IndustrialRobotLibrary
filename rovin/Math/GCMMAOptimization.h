#pragma once

#include "Function.h"
#include <rovin\Utils\Diagnostic.h>

namespace rovin
{
	class GCMMAOptimization
	{
	private:
		int _xN; ///< number of parameters
		int _ineqN; ///< number of inequality constrain

		FunctionPtr _objectFunc;
		FunctionPtr _ineqConstraint;

	public:
		VectorX resultX;
		Real resultFunc;

	public:
		GCMMAOptimization() : _xN(-1), _ineqN(-1),
			_objectFunc(FunctionPtr(new EmptyFunction())), _ineqConstraint(FunctionPtr(new EmptyFunction())) {}
		GCMMAOptimization(int xN, int ineqN) : _xN(xN), _ineqN(ineqN),
			_objectFunc(FunctionPtr(new EmptyFunction())), _ineqConstraint(FunctionPtr(new EmptyFunction())) {
			initialize(xN, ineqN);
		}

		void setXN(int xN) { _xN = xN; }
		void setIneqN(int ineqN) { _ineqN = ineqN; }

		void setObjectiveFunction(const FunctionPtr& objectFunc) { _objectFunc = objectFunc; }
		void setInequalityConstraint(const FunctionPtr& ineqConstraint) { _ineqConstraint = ineqConstraint; }

		FunctionPtr& getObjectiveFunction() { return _objectFunc; }
		FunctionPtr& getInequalityConstraint() { return _ineqConstraint; }


	public:
		void initialize(int xN, int ineqN);
		void setParameters(const Real& albefa, const Real& move0, const Real& asyinit, const Real& asydecr, const Real& asyincr);
		void setCoefficients(const Real& a0, const VectorX& ai, const VectorX& ci, const VectorX& di);
		void setMinMax(const VectorX& minX, const VectorX& maxX);

		void solve(const VectorX& initialX);

		void solveSubProblem(const VectorX& p0, const MatrixX& pi, const VectorX& q0, const MatrixX& qi, const VectorX& ri, /* output */ VectorX& xout);
		void calcLowUpp(int iter, const VectorX& xk, const VectorX& xkm1, const VectorX& xkm2); ///< l^{k}, u^{k}
		void calcAlphaBeta(const VectorX& xk); ///< alpha^{k}, beta^{k}
		void calcPlusMinusMatrix(const MatrixX& mat, /* output */ MatrixX& matp, MatrixX& matm);
		void calcInitialRho(const MatrixX& df0dx, const MatrixX& dfidx, /* output */ Real& rho0, VectorX& rhoi);
		void calcPQR(const Real& rho0, const VectorX& rhoi, const MatrixX& df0dxp, const MatrixX& df0dxm, const MatrixX& dfidxp, const MatrixX& dfidxm,
			const VectorX& xk, const VectorX& f0val, const VectorX& fival, /* output */ VectorX& p0, MatrixX& pi, VectorX& q0, MatrixX& qi, Real& r0, VectorX& ri);
		void initializeSubProb(void);
		void separateFromW(void);
		void separateFromW(const VectorX& w, /* output */ VectorX& x, VectorX& y, Real& z, VectorX& lam, VectorX& xsi, VectorX& eta,
			VectorX& mu, Real& zet, VectorX& s);
		void combineToW(const VectorX& x, const VectorX& y, const Real& z, const VectorX& lam, const VectorX& xsi, const VectorX& eta,
			const VectorX& mu, const Real& zet, const VectorX& s, /* output */ VectorX& w);

		// modification ing
		void calcGradientW(const VectorX& p0, const MatrixX& pi, const VectorX& q0, const MatrixX& qi, const VectorX& bi, /* output */ VectorX& delw);
		void calcGradientW_tmp(const VectorX& p0, const MatrixX& pi, const VectorX& q0, const MatrixX& qi, const VectorX& bi, /* output */ VectorX& delw);

		void calcpqlam(const VectorX& pq0, const MatrixX& pqi, const VectorX& lam, /* output */ VectorX& pqlam);
		void calcdpsidx(const VectorX& plam, const VectorX& qlam, const VectorX& x, /* output */ VectorX& dpsidx);
		void calcgi(const MatrixX& pi, const MatrixX& qi, const VectorX& x, /* output */ VectorX& gival);
		Real calcStepLength(const VectorX& delw, const VectorX& p0, const MatrixX& pi, const VectorX& q0, const MatrixX& qi, const VectorX& bi);


		// modification ing
		Real calcNormResidual(const VectorX& p0, const MatrixX& pi, const VectorX& q0, const MatrixX& qi, const VectorX& bi, const VectorX& delw, Real stepLength, int normCh);
		Real calcNormResidual_tmp(const VectorX& p0, const MatrixX& pi, const VectorX& q0, const MatrixX& qi, const VectorX& bi, const VectorX& delw, Real stepLength, int normCh);



		bool testILSuccess(const VectorX& p0, const MatrixX& pi, const VectorX& q0, const MatrixX& qi, const Real& r0, const VectorX& ri, const VectorX& testx,
			/* output */ VectorX& f0valknu, VectorX& fivalknu, VectorX& f0tvalknu, VectorX& fitvalknu);
		void calcf0tilde(const VectorX& p0, const VectorX& q0, const Real& r0, const VectorX& x, /* output */ VectorX& f0tval);
		void calcfitilde(const MatrixX& pi, const MatrixX& qi, const VectorX& ri, const VectorX& x, /* output */ VectorX& fitval);
		void updateRho0i(const VectorX& xknu, const VectorX& xk, const VectorX& f0valknu, const VectorX& fivalknu, const VectorX& f0tvalknu, const VectorX& fitvalknu, /* output */ Real& rho0, VectorX& rhoi);
	
	
		// steepest gradient method
		//void xfunctionOflam(const VectorX& p0, const MatrixX& pi, const VectorX& q0, const MatrixX& qi,const VectorX& lam, VectorX& x);
		//void yfunctionOflam(const VectorX& lam, VectorX& y);
	
	public:
		// user setting parameters
		Real _ALBEFA;
		Real _MOVE0;
		Real _ASYINIT;
		Real _ASYDECR;
		Real _ASYINCR;

		// coefficients for optimization problem
		Real _a0;
		VectorX _ai;
		VectorX _ci;
		VectorX _di;

		// lower/upper bound for x
		VectorX _minX;
		VectorX _maxX;

		// terminate conditions
		Real _tolX;
		Real _tolFunc;
		int _maxIterOL;
		int _maxIterIL;

	public:
		// variables for outer loop: updated in every outer loop
		VectorX _ollow;
		VectorX _ollowm1;
		VectorX _olupp;
		VectorX _oluppm1;
		VectorX _olalpha;
		VectorX _olbeta;

	public:
		// variables for subproblem
		VectorX _subx;		// _xN
		VectorX _suby;		// _ineqN
		Real	_subz;
		VectorX _sublam;	// _ineqN
		VectorX _subxsi;	// _xN
		VectorX _subeta;	// _xN
		VectorX _submu;		// _ineqN
		Real	_subzet;
		VectorX _subs;		// _ineqN

		VectorX _subw;
		Real	_subeps;
		int		_subDimW;

	public:
		// tmp variable
		VectorX tmpsubx;
		VectorX tmpsublam;
		VectorX tmpplam;
		VectorX tmpqlam;
		VectorX tmpdpsidx;
		VectorX tmpgival;

		MatrixX G;
		MatrixX Gt;
		VectorX xadi;
		VectorX bxdi;
		VectorX ydi;
		VectorX ldi;

		MatrixX Dx;
		MatrixX iDx;
		MatrixX Dly;
		MatrixX iDly;
		VectorX iDy;

		VectorX dxt;
		VectorX dyt;
		Real dzt;
		VectorX dlt;
		VectorX dlyt;

		MatrixX DlyGDxGt; // Dly + GiDxGt , _ineqN by _ineqN
		MatrixX iDlyGDxGt; // inverse of (Dly + GiDxGt), _ineqN by _ineqN
		MatrixX Lower; // Cholesky decomposition lower matrix of DlyGDxGt, _ineqN by _ineqN
		MatrixX iLower; // inverse of lower matrix, _ineqN by _ineqN
		MatrixX iLowert; // transpose of inverse of lower matrix, _ineqN by _ineqN

		MatrixX DxGtiDlyG; // Dx + GtiDlyG, _xN by _xN
		MatrixX iDxGtiDlyG; // inverse of (Dx + GtiDlyG), _xN by _xN
		MatrixX Lower2; // Cholesky decomposition lower matrix of DxGtiDlyG, _xN by _xN
		MatrixX iLower2; // inverse of lower2 matrix, _xN by _xN
		MatrixX iLowert2; // transpose of inverse of lower matrix, _xN by _xN 

		MatrixX MatSizeineqNbyxN;
		MatrixX MatSizexNbyineqN;

		VectorX delx;
		VectorX dellam;

		VectorX resvec;
	};
}