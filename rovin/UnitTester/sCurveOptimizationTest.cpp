#include <iostream>
#include <conio.h>
#include <rovin\Math\Function.h>
#include <memory>
#include <time.h>

#include "rovin\Optimizer\NonlinearOptimization.h"
#include "rovin\Optimizer\UnconstraintedOptimization.h"
#include "rovin\Math\GCMMAOptimization.h"

using namespace std;
using namespace rovin;
using namespace irLib::Opt;

Real makeRandLU(Real lower, Real upper);

//////////////////////// case 1 ////////////////////////
class Step2PosTriZeroNegTri : public Function // unknown : ap1, ap2
{
public:
	Real tsyn, a0, vt, v0, pt, p0, jmax, amax, vmax;

public:
	Step2PosTriZeroNegTri(const Real _tsyn, const Real _a0, const Real _vt, const Real _v0, const Real _pt, const Real _p0, const Real _jmax, const Real _amax, const Real _vmax)
		: tsyn(_tsyn), a0(_a0), vt(_vt), v0(_v0), pt(_pt), p0(_p0), jmax(_jmax), amax(_amax), vmax(_vmax) {}

	VectorX func(const VectorX& x) const
	{
		VectorX val(2);
		Real ap1 = x(0), ap2 = x(1);
		val(0) = -(a0*a0 - 2 * ap1*ap1 + 2 * ap2*ap2 - 2 * jmax*v0 + 2 * jmax*vt) / (2 * jmax);
		val(1) = (6 * a0*ap1*ap1 + 6 * jmax*jmax * p0 - 6 * jmax*jmax * pt - a0*a0*a0 - 6 * ap1*ap1*ap1 +
			6 * ap2*ap2*ap2 + 6 * jmax*jmax * tsyn*v0 - 3 * a0*a0 * jmax*tsyn + 6 * ap1*ap1 * jmax*tsyn) / (6 * jmax*jmax);
		return val;
	}

	MatrixX InverseJacobian(const VectorX & x) const
	{
		MatrixX j(2, 2);
		Real ap1 = x(0), ap2 = x(1);
		Real det = (2 * ap1*ap2*(2 * a0 - 3 * ap1 + 3 * ap2 + 2 * jmax*tsyn)) / (jmax *jmax * jmax);

		j(0, 0) = ((3 * ap2*ap2) / (jmax*jmax)) / det;
		j(1, 1) = ((2 * ap1) / jmax) / det;

		j(0, 1) = ((2 * ap2) / jmax) / det;
		j(1, 0) = -((ap1*(2 * a0 - 3 * ap1 + 2 * jmax*tsyn)) / (jmax*jmax)) / det;

		return j;
	}
};

// case 1 objective function
class case1_objective : public Function
{
public:
	Real tsyn, a0, vt, v0, pt, p0, jmax, amax, vmax;

public:
	case1_objective(const Real _tsyn, const Real _a0, const Real _vt, const Real _v0, const Real _pt, const Real _p0, const Real _jmax, const Real _amax, const Real _vmax)
		: tsyn(_tsyn), a0(_a0), vt(_vt), v0(_v0), pt(_pt), p0(_p0), jmax(_jmax), amax(_amax), vmax(_vmax) {}

	VectorX func(const VectorX& x) const
	{
		VectorX fval(1);
		Real ap1 = x(0), ap2 = x(1), t34 = x(2);
		Real f1 = 0, f2 = 0, f3 = 0;
		f1 = -(a0 - 2 * ap1 + 2 * ap2 - jmax*t34 + jmax*tsyn) / jmax;
		f2 = -(a0 * a0 - 2 * ap1 * ap1 + 2 * ap2 * ap2 - 2 * jmax*v0 + 2 * jmax*vt) / (2 * jmax);
		f3 = (6 * a0 * a0 * ap2 - 6 * a0 * a0 * ap1 - 12 * ap1*ap1 * ap2 + 6 * jmax*jmax * p0 - 
			6 * jmax*jmax * pt + 2 * std::pow(a0, 3) + 6 * std::pow(ap1, 3) + 6 * std::pow(ap2, 3) + 6 * jmax*jmax * t34*v0 -
			6 * a0*jmax*v0 + 12 * ap1*jmax*v0 - 12 * ap2*jmax*v0 - 3 * a0 * a0 * jmax*t34 + 
			6 * ap1*ap1 * jmax*t34) / (6 * jmax*jmax);

		fval(0) = f1*f1 + f2*f2 + f3*f3;

		return fval;
	}

	MatrixX Jacobian(const VectorX& x) const
	{
		MatrixX j(1, 3);
		Real ap1 = x(0), ap2 = x(1), t34 = x(2);
		j(0) = ((12 * jmax*v0 - 24 * ap1*ap2 - 6 * (a0 * a0) + 18 * (ap1 * ap1) + 12 * ap1*jmax*t34)*(6 * (a0 * a0) * ap2 - 6 * (a0 * a0) * ap1 - 
			12 * (ap1 * ap1) * ap2 + 6 * (jmax * jmax) * p0 - 6 * (jmax * jmax) * pt + 2 * std::pow(a0, 3) + 6 * std::pow(ap1, 3) + 6 * std::pow(ap2, 3) + 6 * (jmax * jmax) * t34*v0 - 
			6 * a0*jmax*v0 + 12 * ap1*jmax*v0 - 12 * ap2*jmax*v0 - 3 * (a0 * a0) * jmax*t34 + 6 * (ap1 * ap1) * jmax*t34)) / (18 * std::pow(jmax, 4)) - 
			(4 * a0 - 8 * ap1 + 8 * ap2 - 4 * jmax*t34 + 4 * jmax*tsyn) / (jmax * jmax) - (2 * ap1*((a0 * a0) - 2 * (ap1 * ap1) + 2 * (ap2 * ap2) - 2 * jmax*v0 + 2 * jmax*vt)) / (jmax * jmax);
		j(1) = (4 * a0 - 8 * ap1 + 8 * ap2 - 4 * jmax*t34 + 4 * jmax*tsyn) / (jmax * jmax) + (2 * ap2*((a0 * a0) - 2 * (ap1 * ap1) + 
			2 * (ap2 * ap2) - 2 * jmax*v0 + 2 * jmax*vt)) / (jmax * jmax) - ((-6 * (a0 * a0) + 12 * (ap1 * ap1) - 18 * (ap2 * ap2) + 
				12 * jmax*v0)*(6 * (a0 * a0) * ap2 - 6 * (a0 * a0) * ap1 - 12 * (ap1 * ap1) * ap2 + 6 * (jmax * jmax) * p0 - 6 * (jmax * jmax) * pt + 
					2 * std::pow(a0, 3) + 6 * std::pow(ap1, 3) + 6 * std::pow(ap2, 3) + 6 * (jmax * jmax) * t34*v0 - 6 * a0*jmax*v0 + 12 * ap1*jmax*v0 - 12 * ap2*jmax*v0 - 
					3 * (a0 * a0) * jmax*t34 + 6 * (ap1 * ap1) * jmax*t34)) / (18 * std::pow(jmax, 4));
		j(2) = ((-(a0 * a0) + 2 * (ap1 * ap1) + 2 * jmax*v0)*(6 * (a0 * a0) * ap2 - 6 * (a0 * a0) * ap1 - 12 * (ap1 * ap1) * ap2 + 
			6 * (jmax * jmax) * p0 - 6 * (jmax * jmax) * pt + 2 * std::pow(a0, 3) + 6 * std::pow(ap1, 3) + 6 * std::pow(ap2, 3) + 6 * (jmax * jmax) * t34*v0 - 
			6 * a0*jmax*v0 + 12 * ap1*jmax*v0 - 12 * ap2*jmax*v0 - 3 * (a0 * a0) * jmax*t34 + 6 * (ap1 * ap1) * jmax*t34)) / (6 * std::pow(jmax, 3)) -
			(2 * (a0 - 2 * ap1 + 2 * ap2 - jmax*t34 + jmax*tsyn)) / jmax;

		return j;
	}


};

// case 1 inequality constraint function
class case1_constraint : public Function
{
public:
	Real tsyn, a0, vt, v0, pt, p0, jmax, amax, vmax;

public:
	case1_constraint(const Real _tsyn, const Real _a0, const Real _vt, const Real _v0, const Real _pt, const Real _p0, const Real _jmax, const Real _amax, const Real _vmax)
		: tsyn(_tsyn), a0(_a0), vt(_vt), v0(_v0), pt(_pt), p0(_p0), jmax(_jmax), amax(_amax), vmax(_vmax) {}

public:
	VectorX func(const VectorX& x) const
	{
		VectorX fval(5);
		Real ap1 = x(0), ap2 = x(1), t34 = x(2);
		fval(0) = -ap1;
		fval(1) = ap1 - amax;
		fval(2) = ap2;
		fval(3) = -ap2 - amax;
		fval(4) = -t34;
		return fval;
	}

	MatrixX Jacobian(const VectorX& x) const
	{
		MatrixX j(3, 5);
		j.setZero();
		j(0, 0) = -1;
		j(0, 1) = 1;
		j(1, 2) = 1;
		j(1, 3) = -1;
		j(2, 4) = -1;
		return j;
		//VectorX j(5, 3);
		//j.setZero();
		//j(0, 0) = -1;
		//j(1, 0) = 1;
		//j(2, 1) = 1;
		//j(3, 1) = -1;
		//j(4, 2) = -1;
		//return j;
	}
};

// case 1 inequality constraint function
class case1_constraint_GCMMA : public Function
{
public:
	case1_constraint_GCMMA(){}
	VectorX func(const VectorX& x) const
	{
		VectorX fval(1);
		fval(0) = -1;
		return fval;
	}

	MatrixX Jacobian(const VectorX& x) const
	{
		MatrixX j(1, 3);
		j(0, 0) = 0; j(0, 1) = 0; j(0, 2) = 0;
		return j;
	}
};

int main()
{
	// initial value
	VectorX vec2(2), vec3(3);
	srand(time(NULL));
	Real value2 = rand();
	NewtonRaphson solver(3, 3);
	Real time;

	// PosTriZeroNegTri Case Experiment --> numerical
	cout << "------- case 1 : PosTriZeroNegTri -------" << endl;

	// newton raphon method
	int numOfExp = 10;
	int errorcnt = 0;
	bool r = true;
	for (int k = 0; k < numOfExp; k++)
	{
		r = true;
		std::shared_ptr<Step2PosTriZeroNegTri> obj1 = std::shared_ptr<Step2PosTriZeroNegTri>(new Step2PosTriZeroNegTri(3, 5, 6.25, 0, 19.7917, 0, 10, 100, 100));
		vec2(0) = makeRandLU(0, 10); vec2(1) = makeRandLU(-10, 0);
		solver.setxN(2); solver.setfN(2); solver.setfunction(obj1);
		solver.solve(vec2);
		if (std::pow((10 - solver.getResultX()(0)), 2) > 1E-5)
		{
			errorcnt++;
			r = false;
		}
		cout << "initial ap1 : " << vec2(0) << ", initial ap2 : " << vec2(1) << "/ Result : " << solver.getResultX()(0) << ", " << solver.getResultX()(1) << ", " << r << endl;
	}

	Real success = (Real)(numOfExp - errorcnt) / (Real)(numOfExp)* 100;
	cout << "Newton raphson success per : " << success << "%" << endl << endl;

	// constraint optimization method
	// nlopt & GCMMA
	std::shared_ptr<case1_objective> obj = std::shared_ptr<case1_objective>(new case1_objective(3, 5, 6.25, 0, 19.7917, 0, 10, 30, 30));
	std::shared_ptr<case1_constraint> con = std::shared_ptr<case1_constraint>(new case1_constraint(3, 5, 6.25, 0, 19.7917, 0, 10, 30, 30));
	NonlinearOptimization nonSolver(3, 0, 5);
	nonSolver.setObjectiveFunction(obj);
	nonSolver.setInequalityConstraint(con);
	//VectorX init(3); init(0) = 20; init(1) = -20; init(2) = 5;
	VectorX init(3); init(0) = 12; init(1) = -6; init(2) = 1;

	cout << obj->func(init) << endl;
	cout << obj->Jacobian(init) << endl;

	time = clock();
	nonSolver.solve(init);
	cout << "computation time : " << clock() - time << endl << endl;
	cout << nonSolver.resultX << endl;
	cout << nonSolver.resultFunc << endl << endl;

	//cout << "---------------GCMMA---------------" << endl;
	//
	//GCMMAOptimization* opt = new GCMMA_PDIPM(3, 1);
	//VectorX  minX(3), maxX(3);
	//minX(0) = 0; minX(1) = -30; minX(2) = 0;
	//maxX(0) = 30; maxX(1) = 0; maxX(2) = 15;
	//
	//std::shared_ptr<case1_constraint_GCMMA> con_GCMMA = std::shared_ptr<case1_constraint_GCMMA>(new case1_constraint_GCMMA());

	//opt->setObjectiveFunction(obj);
	//opt->setInequalityConstraint(con_GCMMA);
	//opt->setMinMax(minX, maxX);
	//opt->setTolX(1E-7);
	//opt->setTolFunc(1E-7);
	////opt->setMaxIterOL(3000);
	//GCMMAReturnFlag retFlag;
	//time = clock();
	//retFlag = opt->solve(init);
	//cout << "computation time : " << clock() - time << endl << endl;
	//cout << opt->getResultX() << endl;
	//cout << opt->getResultFunc() << endl << endl;


	cout << "Program Complete" << endl;
	_getch();
	return 0;
}

Real makeRandLU(Real lower, Real upper)
{
	return lower + (upper - lower)*((double)rand() / 32767.0);
}