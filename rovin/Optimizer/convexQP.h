#pragma once

#include "rovin\Math\Common.h"
#include "rovin\Math\Constant.h"
#include "rovin\Utils\Diagnostic.h"
#include <list>

namespace rovin
{
	namespace irOpt
	{
		class convexQP;


		class convexQP
		{
		public:
			convexQP(){}
			convexQP(const MatrixX& G, const VectorX& c, const MatrixX& A, const VectorX& b) :_G(G), _c(c), _A(A), _b(b) { _xN = c.rows(); _ineqN = b.rows(); }
			convexQP(int xN, int ineqN, const MatrixX& G, const VectorX& c, const MatrixX& A, const VectorX& b) :_xN(xN), _ineqN(ineqN), _G(G), _c(c), _A(A), _b(b) {}
			~convexQP() { _idxWS.clear(); _idxWSdual.clear(); }

			void setDimension(int xN, int ineqN) { _xN = xN; _ineqN = ineqN; }
			void setProblem(const MatrixX& G, const VectorX& c, const MatrixX& A, const VectorX& b) { _G = G; _c = c; _A = A; _b = b; _xN = c.rows(); _ineqN = b.rows(); }
			void checkDimension(void)
			{
				LOGIF(_G.rows() == _xN,		"dimension error at _G.rows() or _xN");
				LOGIF(_G.cols() == _xN,		"dimension error at _G.cols() or _xN");
				LOGIF(_c.rows() == _xN,		"dimension error at _c or _xN");
				LOGIF(_A.rows() == _ineqN,	"dimension error at _A.rows() or _ineqN");
				LOGIF(_A.cols() == _xN,		"dimension error at _A.cols() or _xN");
				LOGIF(_b.rows() == _ineqN,	"dimension error at _b or _ineqN")
			}

		public:

			///> solve " minimize 0.5 * x^T G x + c^T x   subject to A*x >= b
			bool solveInequalityConstrainedQP(void);
			void initialize(void);
			bool checkSuccess(void);
			bool checkFailure(void);
			void selectIndexR(void);
			void calcPandV(void);
			void calcTandIndex(void);
			void updateMu(const Real& t);
			void appendConToWS(int idx);
			void removeConFromWS(int idx);

			static Real solveEqualityConstrainedQP(const MatrixX& G, const VectorX& c, const MatrixX& A, const VectorX& b, VectorX& xsol);
		
		private:
			MatrixX _G;
			VectorX _c;
			MatrixX _A;
			VectorX _b;
			VectorX _ineqCon;
			int _xN;
			int _ineqN;

		public:
			VectorX _resultX;	///> dim: _xN
			VectorX _resultmu;	///> dim: _ineqN
			Real _resultFunc;	///> result value of objective function

		private:
			VectorX _x;
			VectorX _mu;

			int _noWS;					///> number of constraints in working set
			int _noWSdual;				///> number of constraints in dual working set
			std::list<int> _idxWS;		///> indices of working set
			std::list<int> _idxWSdual;	///> indices of dual working set
			std::list<int>::iterator _iterList;

			int _idxr;	///> index of inequality that has least negative value in idxWSdual

			VectorX _p; ///> dim: _xN;
			VectorX _v; ///> dim varies in every loop, depends on noWS
			MatrixX _tmpA; VectorX _tmpx, _tmpb; ///> variables for calculating p and v
			
			Real _t1, _t2, _t3; ///> candidate step length
			int _idx1; ///> candidate index of constraint to be removed
		};
	}
}
