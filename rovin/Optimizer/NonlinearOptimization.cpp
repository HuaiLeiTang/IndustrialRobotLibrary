#include "NonlinearOptimization.h"

#include <nlopt.hpp>
#include <rovin\Math\Common.h>

using namespace std;
using namespace nlopt;

namespace rovin
{
	double objective(unsigned n, const double* x, double* grad, void *f_data)
	{
		NonlinearOptimization* ptr = reinterpret_cast<NonlinearOptimization*>(f_data);
		VectorX xValue(n);
		for (unsigned int i = 0; i < n; i++)
		{
			xValue(i) = x[i];
		}
		if (grad)
		{
			MatrixX jacobian = (*ptr->getObjectiveFunction()).Jacobian(xValue);
			for (unsigned int j = 0; j < n; j++)
			{
				grad[j] = jacobian(0, j);
			}
		}
		return (*ptr->getObjectiveFunction())(xValue)(0);
	}

	void constraint(unsigned m, double* result, unsigned n, const double* x, double* grad, void* f_data)
	{
		pair<bool, NonlinearOptimization*>* data = reinterpret_cast<pair<bool, NonlinearOptimization*>*>(f_data);
		FunctionPtr constraintFunction;
		if (!data->first)
		{
			constraintFunction = data->second->getEqualityConstraint();
		}
		else
		{
			constraintFunction = data->second->getInequalityConstraint();
		}
		VectorX xValue(n);
		for (unsigned int i = 0; i < n; i++)
		{
			xValue(i) = x[i];
		}
		if (grad)
		{
			MatrixX jacobian = (*constraintFunction).Jacobian(xValue);
			for (unsigned int i = 0; i < m; i++)
			{
				for (unsigned int j = 0; j < n; j++)
				{
					grad[i*n + j] = jacobian(i, j);
				}
			}
		}
		VectorX fval = (*constraintFunction).func(xValue);
		for (unsigned int i = 0; i < m; i++)
		{
			result[i] = fval(i);
		}
	}

	void NonlinearOptimization::solve(const VectorX& initialX)
	{
		if (_xN == -1)		_xN = initialX.size();
		if (_eqN == -1)		_eqN = _eqConstraint->func(initialX).size();
		if (_ineqN == -1)	_ineqN = _ineqConstraint->func(initialX).size();

		opt optimizer;

		if (_eqN == 0 && _ineqN == 0)
		{
			optimizer = opt(nlopt::LD_LBFGS, _xN);
			//optimizer = opt(nlopt::LD_TNEWTON_PRECOND, _xN);
			//optimizer = opt(nlopt::LD_VAR1, _xN);
		}
		else if(_eqN != 0) optimizer = opt(nlopt::LD_SLSQP, _xN);
		else optimizer = opt(nlopt::LD_MMA, _xN);

		optimizer.set_min_objective(objective, this);
		//optimizer.set_min_objective(objective, NULL);
		if (_eqN != 0)
		{
			optimizer.add_equality_mconstraint(constraint, new pair<bool, NonlinearOptimization*>(false, this), vector<Real>(_eqN, 1e-7));
		}
		if (_ineqN != 0)
		{
			optimizer.add_inequality_mconstraint(constraint, new pair<bool, NonlinearOptimization*>(true, this), vector<Real>(_ineqN, 1e-7));
		}

		optimizer.set_xtol_rel(1e-4);
		//optimizer.set_xtol_rel(1e-8);
		optimizer.set_ftol_rel(1e-4);
		//optimizer.set_ftol_rel(1e-8);
		optimizer.set_maxeval(500);
		//optimizer.set_maxeval(3000);

		std::vector<Real> x(_xN);
		for (int i = 0; i < _xN; i++)
		{
			x[i] = initialX(i);
		}

		/////////////////////////////
		//optimizer.optimize(x, resultFunc);

		//resultX = VectorX::Zero(_xN);
		//for (int i = 0; i < _xN; i++)
		//{
		//	resultX(i) = x[i];
		//}
		/////////////////////////////

		VectorX minX = initialX;
		Real min = RealMax;
		for (unsigned int i = 0; i < 10; i++) // 3 -> 10
		{
			Real past = resultFunc;
			bool flag = false;
			try
			{
				optimizer.optimize(x, resultFunc);
			}
			catch (exception&)
			{
				//cout << "[OPTIMIZATION] FORCE STOP" << endl;
			}

			VectorX x_tmp(_xN);
			for (int i = 0; i < _xN; i++)
			{
				x_tmp(i) = x[i];
			}
			resultFunc = _objectFunc->func(x_tmp)(0);

			if (min > resultFunc && RealLessEqual(_ineqConstraint->func(x_tmp), 1e-3))
			{
				minX = x_tmp;
				min = resultFunc;
			}

			if (RealLess(Abs(past - resultFunc), resultFunc / (Real)1.0e+3))
			{
				if (RealEqual(min, RealMax))
				{
					minX = x_tmp;
					min = resultFunc;
				}
				break;
			}
		}

		if (!RealLessEqual(_ineqConstraint->func(minX), 1e-3))
		{
			min = RealMax;
		}

		//cout << "[OPTIMIZATION] Computation Time : " << clock() - startTime << "ms" << endl;

		/////////////////////////////// 만약에 정답이 없는 경우 처리 ///////////////////////////////
		resultX = minX;
		resultFunc = min;


	}
}