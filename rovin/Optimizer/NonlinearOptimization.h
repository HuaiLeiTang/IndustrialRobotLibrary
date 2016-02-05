#pragma once

#include <rovin\Math\Function.h>

namespace rovin
{
	class NonlinearOptimization
	{
	public:
		NonlinearOptimization() : _xN(-1), _eqN(-1), _ineqN(-1), 
			_objectFunc(FunctionPtr(new EmptyFunction())), _eqConstraint(FunctionPtr(new EmptyFunction())), _ineqConstraint(FunctionPtr(new EmptyFunction())) {}
		NonlinearOptimization(int xN, int eqN, int ineqN) : _xN(xN), _eqN(eqN), _ineqN(ineqN),
			_objectFunc(FunctionPtr(new EmptyFunction())), _eqConstraint(FunctionPtr(new EmptyFunction())), _ineqConstraint(FunctionPtr(new EmptyFunction())) {}

		void setXN(int xN) { _xN = xN; }
		void setEqN(int eqN) { _eqN = eqN; }
		void setIneqN(int ineqN) { _ineqN = ineqN; }

		void setObjectiveFunction(const FunctionPtr& objectFunc) { _objectFunc = objectFunc; }
		void setEqualityConstraint(const FunctionPtr& eqConstraint) { _eqConstraint = eqConstraint; }
		void setInequalityConstraint(const FunctionPtr& ineqConstraint) { _ineqConstraint = ineqConstraint; }

		FunctionPtr& getObjectiveFunction() { return _objectFunc; }
		FunctionPtr& getEqualityConstraint() { return _eqConstraint; }
		FunctionPtr& getInequalityConstraint() { return _ineqConstraint; }

		void solve(const VectorX& initialX);

		VectorX resultX;
		Real resultFunc;

	private:
		int _xN, _eqN, _ineqN;

		FunctionPtr _objectFunc;
		FunctionPtr _eqConstraint;
		FunctionPtr _ineqConstraint;
	};
}