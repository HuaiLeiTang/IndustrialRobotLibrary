#pragma once

#include "Function.h"
#include <rovin\Utils\Diagnostic.h>


//#define STRATEGY_01
//#define STRATEGY_02


namespace rovin
{
	class GCMMAOptimization
	{
	public:
		int _xN; ///< number of parameters
		int _ineqN; ///< number of inequality constrain
		FunctionPtr _objectFunc;
		FunctionPtr _ineqConstraint;

		// coefficients for optimization problem
		Real _a0;
		VectorX _ai;
		VectorX _ci;
		VectorX _di;

		// variables for outer loop: updated in every outer loop
		VectorX _ollow;
		VectorX _ollowm1;
		VectorX _olupp;
		VectorX _oluppm1;
		VectorX _olalpha;
		VectorX _olbeta;

		// variables for inner loop: updated in every inner loop
		VectorX _ilp0;
		VectorX _ilq0;
		MatrixX _ilpi;
		MatrixX _ilqi;
		Real _ilr0;
		Real _ilrho0;
		VectorX _ilri;
		VectorX _ilrhoi;

#ifdef STRATEGY_01
		VectorX _olsigma;
		VectorX _ols;
		VectorX _olt;
		Real _oleta;
		VectorX _olb;
#endif

	public:
		GCMMAOptimization() : _xN(-1), _ineqN(-1),
			_objectFunc(FunctionPtr(new EmptyFunction())), _ineqConstraint(FunctionPtr(new EmptyFunction())) {}
		GCMMAOptimization(int xN, int ineqN) : _xN(xN), _ineqN(ineqN),
			_objectFunc(FunctionPtr(new EmptyFunction())), _ineqConstraint(FunctionPtr(new EmptyFunction())) {
			initialize(xN, ineqN);
		}

		

		// get functions
		const FunctionPtr& getObjectiveFunction() const { return _objectFunc; }
		const FunctionPtr& getInequalityConstraint() const { return _ineqConstraint; }
		FunctionPtr& getObjectiveFunction() { return _objectFunc; }
		FunctionPtr& getInequalityConstraint() { return _ineqConstraint; }
		const VectorX& getResultX() const { return resultX; }
		const Real getResultFunc() const { return resultFunc; }

		// set functions
		void setXN(const int xN) { _xN = xN; }
		void setIneqN(const int ineqN) { _ineqN = ineqN; }
		void setObjectiveFunction(const FunctionPtr& objectFunc) { _objectFunc = objectFunc; }
		void setInequalityConstraint(const FunctionPtr& ineqConstraint) { _ineqConstraint = ineqConstraint; }
		void setParameters(const Real& albefa, const Real& move0, const Real& asyinit, const Real& asydecr, const Real& asyincr);
		void setCoefficients(const Real& a0, const VectorX& ai, const VectorX& ci, const VectorX& di);
		void setMinMax(const VectorX& minX, const VectorX& maxX);


		void solve(const VectorX& initialX);

		virtual ~GCMMAOptimization() {}

	public:
		

#ifdef STRATEGY_01
		void calcSigma(int iter, const VectorX& xk, const VectorX& xkm1, const VectorX& xkm2);
		void calcInitialRho_st01(int iter, const VectorX& xk, const VectorX& xkm1, const MatrixX& df0dx, const MatrixX& df0dxm1,
			const MatrixX& dfidx, const MatrixX& dfidxm1);
#endif

		void calcInitialRho(const MatrixX& df0dx, const MatrixX& dfidx);
		void calcPQR(const MatrixX& df0dxp, const MatrixX& df0dxm, const MatrixX& dfidxp, const MatrixX& dfidxm, const VectorX& xk, const VectorX& f0val, const VectorX& fival);
		void updateRho0i(const VectorX& xknu, const VectorX& xk, const VectorX& f0valknu, const VectorX& fivalknu, const VectorX& f0tvalknu, const VectorX& fitvalknu);
		bool testILSuccess(const VectorX& testx, /* output */ VectorX& f0valknu, VectorX& fivalknu, VectorX& f0tvalknu, VectorX& fitvalknu);
		void calcf0tilde(const VectorX& x, /* output */ VectorX& f0tval);
		void calcfitilde(const VectorX& x, /* output */ VectorX& fitval);

		virtual void solveSubProblem(/* output */ VectorX& xout) = 0;

		// memory allocation for variables of outer loop/inner loop/sub problem
		void allocOLvar(void);
		void allocILvar(void);
		virtual void allocSUBvar(void) = 0;

		void initialize(int xN, int ineqN);

	private:
		void calcLowUpp(int iter, const VectorX& xk, const VectorX& xkm1, const VectorX& xkm2); ///< l^{k}, u^{k}
		void calcAlphaBeta(const VectorX& xk); ///< alpha^{k}, beta^{k}
		void calcPlusMinusMatrix(const MatrixX& mat, /* output */ MatrixX& matp, MatrixX& matm);

	private:
		VectorX resultX;
		Real resultFunc;

		// user setting parameters
		Real _ALBEFA;
		Real _MOVE0;
		Real _ASYINIT;
		Real _ASYDECR;
		Real _ASYINCR;

		// lower/upper bound for x
		VectorX _minX;
		VectorX _maxX;

		// terminate conditions
		Real _tolX;
		Real _tolFunc;
		int _maxIterOL;
		int _maxIterIL;
	
	};

	class GCMMA_TRM : public GCMMAOptimization
	{
		// solve GCMMA subproblem by using 'Trust-Region Method'
	public:
		GCMMA_TRM() : GCMMAOptimization() { setParametersTR(0.4, 0.6, 0.2, 0.4, 1.5); }
		GCMMA_TRM(int xN, int ineqN) : GCMMAOptimization(xN, ineqN) { setParametersTR(0.3, 0.7, 0.5, 0.7, 1.2); }

	public:
		void solveSubProblem(/* output */ VectorX& xout);
		void allocSUBvar(void);

		void setParametersTR(const Real& v, const Real& w, const Real& gam0, const Real& gam1, const Real& gam2);
		void initializeSubProb(void);

		void calcx(const VectorX& lam, /* output */ VectorX& subx);
		void calcy(const VectorX& lam, /* output */ VectorX& suby);
		void calcW(const VectorX& lam, /* output */ Real& W);
		void calcdW(const VectorX& lam, /* output */ VectorX& dW);
		void calcm(const VectorX& lam, /* output */ Real& m);
		void calceta(void);
		void calclamhat(void);

	public:
		// parameters for trust-region
		Real _TR_v;
		Real _TR_w;
		Real _TR_gamma0;
		Real _TR_gamma1;
		Real _TR_gamma2;

		VectorX _sublam, _sublamm1, _sublamhat;
		VectorX _subdW, _subdWm1;
		Real _subW, _subWhat;
		Real _subm, _submhat;

		VectorX _subs; // s = lam - lamm1;
		VectorX _subt; // t = dw - dwm1;
		Real _subeta;
		Real _radius;
		Real _subtheta;

		VectorX _subx;
		VectorX _suby;
	};

	class GCMMA_PDIPM : public GCMMAOptimization
	{
		// solve GCMMA subproblem by using 'Primal-Dual Interior-Point Method'
	public:
		GCMMA_PDIPM() : GCMMAOptimization() {}
		GCMMA_PDIPM(int xN, int ineqN) : GCMMAOptimization(xN, ineqN) {}

	public:
		void solveSubProblem(/* output */ VectorX& xout);
		void allocSUBvar(void);

		void initializeSubProb(void);
		void separateFromW(void);
		void separateFromW(const VectorX& w, /* output */ VectorX& x, VectorX& y, Real& z, VectorX& lam, VectorX& xsi, VectorX& eta,
			VectorX& mu, Real& zet, VectorX& s);
		void combineToW(const VectorX& x, const VectorX& y, const Real& z, const VectorX& lam, const VectorX& xsi, const VectorX& eta,
			const VectorX& mu, const Real& zet, const VectorX& s, /* output */ VectorX& w);
		// modification ing
		//      void calcGradientW(const VectorX& p0, const MatrixX& pi, const VectorX& q0, const MatrixX& qi, const VectorX& bi, /* output */ VectorX& delw);
		void calcGradientW_tmp(/* output */ VectorX& delw);
		void calcpqlam(const VectorX& pq0, const MatrixX& pqi, const VectorX& lam, /* output */ VectorX& pqlam);
		void calcdpsidx(const VectorX& plam, const VectorX& qlam, const VectorX& x, /* output */ VectorX& dpsidx);
		void calcgi(const VectorX& x, /* output */ VectorX& gival);
		Real calcStepLength(const VectorX& delw);
		// modification ing
		//      Real calcNormResidual(const VectorX& p0, const MatrixX& pi, const VectorX& q0, const MatrixX& qi, const VectorX& bi, const VectorX& delw, Real stepLength, int normCh);
		Real calcNormResidual_tmp(const VectorX& delw, Real stepLength, int normCh);

	public:
		// variables for subproblem
		VectorX _subx;      // _xN
		VectorX _suby;      // _ineqN
		Real   _subz;
		VectorX _sublam;   // _ineqN
		VectorX _subxsi;   // _xN
		VectorX _subeta;   // _xN
		VectorX _submu;      // _ineqN
		Real   _subzet;
		VectorX _subs;      // _ineqN

		VectorX _subw;
		Real   _subeps;
		int      _subDimW;

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