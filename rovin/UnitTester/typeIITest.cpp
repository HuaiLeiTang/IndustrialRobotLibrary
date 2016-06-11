#include <iostream>
#include <conio.h>
#include <time.h>
#include "rovin\Math\Function.h"
#include "rovin\Optimizer\NonlinearOptimization.h"
#include "rovin\Optimizer\UnconstraintedOptimization.h"

using namespace std;
using namespace rovin;
using namespace irLib::Opt;

class Objective : public Function
{
public:
	Real vt, v0, jmax, a0, pt, p0;
public:
	Objective(const Real _vt, const Real _v0, const Real _jmax, const Real _a0, const Real _pt, const Real _p0) : 
		vt(_vt), v0(_v0), jmax(_jmax), a0(_a0), pt(_pt), p0(_p0){}

	VectorX func(const VectorX& x) const
	{
		Real ap1 = x(0), ap2 = x(1);
		VectorX val(1);
		val(0) = 0;
		Real f1 = 0, f2 = 0;
		//f1 = vt - v0 - (ap1*ap1 / jmax + 0.5*(a0 - 1)*ap1 / jmax - 0.5*a0*a0 / jmax - ap2*ap2 / jmax);
		f1 = vt - v0 - (ap1*ap1 / jmax - 0.5*a0*a0 / jmax - ap2*ap2 / jmax);
		f2 = pt - p0 - (-2 * vt*ap2 / jmax - ap2*ap2*ap2 / (jmax*jmax) + 2*v0*ap1 / jmax - v0*a0 / jmax + a0*a0*a0 / (3*jmax*jmax) + ap1*ap1*ap1 / (jmax*jmax) - a0*a0*ap1 / (jmax*jmax));
		val(0) = f1*f1 + f2*f2;
		return val;
	}

	void result(Real ap1, Real ap2)
	{
		Real f1 = 0, f2 = 0;
		f1 = ap1*ap1 / jmax - 0.5*a0*a0 / jmax - ap2*ap2 / jmax;
		f2 = -2 * vt*ap2 / jmax - ap2*ap2*ap2 / (jmax*jmax) + 2 * v0*ap1 / jmax - v0*a0 / jmax + a0*a0*a0 / (3 * jmax*jmax) + ap1*ap1*ap1 / (jmax*jmax) - a0*a0*ap1 / (jmax*jmax);
		cout << "vt - v0 : " << vt - v0 << endl;
		cout << "pt - p0 : " << pt - p0 << endl;
		cout << "f1 : " << f1 << endl;
		cout << "f2 : " << f2 << endl;
	}
};

class fcn : public Function
{
public:
	Real vt, v0, jmax, a0, pt, p0;
public:
	fcn(const Real _vt, const Real _v0, const Real _jmax, const Real _a0, const Real _pt, const Real _p0) :
		vt(_vt), v0(_v0), jmax(_jmax), a0(_a0), pt(_pt), p0(_p0) {}

	VectorX func(const VectorX& x) const
	{
		Real ap1 = x(0), ap2 = x(1);
		VectorX val(2);
		val(0) = vt - v0 - (ap1*ap1 / jmax - 0.5*a0*a0 / jmax - ap2*ap2 / jmax);
		val(1) = pt - p0 - (-2 * vt*ap2 / jmax - ap2*ap2*ap2 / (jmax*jmax) + 2 * v0*ap1 / jmax - v0*a0 / jmax + a0*a0*a0 / (3 * jmax*jmax) + ap1*ap1*ap1 / (jmax*jmax) - a0*a0*ap1 / (jmax*jmax));
		return val;
	}

	MatrixX InverseJacobian(const VectorX & x) const
	{
		MatrixX j(2, 2);
		Real ap1 = x(0), ap2 = x(1);
		Real det = -(2 * (a0*a0 * ap2 - 3 * ap1*ap1 * ap2 + 3 * ap1*ap2*ap2 + 2 * jmax*vt*ap1 - 2 * jmax*v0*ap2)) / (jmax*jmax*jmax);

		j(0, 0) = ((2 * vt) / jmax + (3 * ap2 * ap2) / (jmax*jmax)) / det;
		j(0, 1) = -((2 * ap2) / jmax) / det;
		j(1, 0) = -(a0 * a0 / (jmax *jmax) - (2 * v0) / jmax - (3 * ap1 * ap1) / (jmax*jmax)) / det;
		j(1, 1) = (-(2 * ap1) / jmax) / det;

		return j;
	}
};

class case1 : public Function // PosTriNegTri ( a >= 0), output ap1, ap2
{
public:
	Real vt, v0, jmax, a0, pt, p0;
public:
	case1(const Real _vt, const Real _v0, const Real _jmax, const Real _a0, const Real _pt, const Real _p0) :
		vt(_vt), v0(_v0), jmax(_jmax), a0(_a0), pt(_pt), p0(_p0) {}

	VectorX func(const VectorX& x) const
	{
		Real ap1 = x(0), ap2 = x(1);
		VectorX val(2);
		val(0) = vt - v0 - (ap1*ap1 / jmax - 0.5*a0*a0 / jmax - ap2*ap2 / jmax);
		val(1) = pt - p0 - (-2 * vt*ap2 / jmax - ap2*ap2*ap2 / (jmax*jmax) + 2 * v0*ap1 / jmax - v0*a0 / jmax + a0*a0*a0 / (3 * jmax*jmax) + ap1*ap1*ap1 / (jmax*jmax) - a0*a0*ap1 / (jmax*jmax));
		return val;
	}

	MatrixX InverseJacobian(const VectorX & x) const
	{
		MatrixX j(2, 2);
		Real ap1 = x(0), ap2 = x(1);
		Real det = -(2 * (a0*a0 * ap2 - 3 * ap1*ap1 * ap2 + 3 * ap1*ap2*ap2 + 2 * jmax*vt*ap1 - 2 * jmax*v0*ap2)) / (jmax*jmax*jmax);

		j(0, 0) = ((2 * vt) / jmax + (3 * ap2 * ap2) / (jmax*jmax)) / det;
		j(0, 1) = -((2 * ap2) / jmax) / det;
		j(1, 0) = -(a0 * a0 / (jmax *jmax) - (2 * v0) / jmax - (3 * ap1 * ap1) / (jmax*jmax)) / det;
		j(1, 1) = (-(2 * ap1) / jmax) / det;

		return j;
	}
};

class case2 : public Function // PosTriZeroNegTri ( a >= 0 ), output : ap1, ap2, e34
{
public:
	Real a0, vt, v0, pt, p0, jmax, amax, vmax;

public:
	case2(const Real _a0, const Real _vt, const Real _v0, const Real _pt, const Real _p0, const Real _jmax, const Real _amax, const Real _vmax) 
		: a0(_a0), vt(_vt), v0(_v0), pt(_pt), p0(_p0), jmax(_jmax), amax(_amax), vmax(_vmax) {}

	VectorX func(const VectorX& x) const
	{
		Real ap1 = x(0), ap2 = x(1), e23 = x(2);
		VectorX val(3);
		val(0) = -(a0 * a0 - 2 * ap1 * ap1 + 2 * ap2 * ap2 - 2 * jmax*v0 + 2 * jmax*vt) / (2 * jmax);
		val(1) = (-a0 * a0 + 2 * ap1 * ap1 + 2 * jmax*v0 - 2 * jmax*vmax) / (2 * jmax);
		val(2) = -(6 * a0 * a0 * ap1 - 3 * a0 * a0 * ap2 + 6 * ap1 * ap1 * ap2 - 6 * jmax * jmax * p0 + 6 * jmax * jmax * pt -
			2 * a0 * a0 * a0 - 6 * ap1 * ap1 * ap1 - 6 * e23 * jmax * jmax * v0 + 6 * a0*jmax*v0 - 12 * ap1*jmax*v0 +
			6 * ap2 * jmax*v0 + 6 * ap2 * jmax*vt + 3 * a0 * a0 * e23 * jmax - 6 * ap1 * ap1 * e23 * jmax) / (6 * jmax * jmax);
		return val;
	}

	MatrixX InverseJacobian(const VectorX & x) const
	{
		MatrixX j(3, 3); j.setZero();
		Real ap1 = x(0), ap2 = x(1), e23 = x(2);
		j(0, 1) = jmax / (2 * ap1);
		j(1, 0) = -jmax / (2 * ap2);
		j(1, 1) = jmax / (2 * ap2);
		j(2, 0) = -(-a0 * a0 + 2 * ap1 * ap1 + 2 * jmax*v0 + 2 * jmax*vt) / (2 * ap2*(-a0 * a0 + 2 * ap1 * ap1 + 2 * jmax*v0));
		j(2, 1) = (2 * a0 * a0 * ap2 - a0 * a0 * ap1 + 4 * ap1*ap2*ap2 - 6 * ap1 *ap1 * ap2 + 2 * ap1 *ap1 * ap1 + 2 * ap1*jmax*v0 - 4 * ap2*jmax*v0 + 2 * ap1*jmax*vt - 4 * ap1*ap2*e23*jmax) / (2 * ap1*ap2*(-a0*a0 + 2 * ap1*ap1 + 2 * jmax*v0));
		j(2, 2) = (2 * jmax) / (-a0 * a0 + 2 * ap1 * ap1 + 2 * jmax*v0);
		return j;
	}

};

class case3 : public Function // PosTrapNegTri ( a >= 0 ), output : ap2, e12
{
public:
	Real a0, vt, v0, pt, p0, jmax, amax, vmax;

public:
	case3(const Real _a0, const Real _vt, const Real _v0, const Real _pt, const Real _p0, const Real _jmax, const Real _amax, const Real _vmax)
		: a0(_a0), vt(_vt), v0(_v0), pt(_pt), p0(_p0), jmax(_jmax), amax(_amax), vmax(_vmax) {}

	VectorX func(const VectorX& x) const
	{
		Real ap2 = x(0), e12 = x(1);
		VectorX val(2);
		val(0) = (-a0 *a0 + 2 * amax * amax + 2 * e12*jmax*amax - 2 * ap2 * ap2 + 2 * jmax*v0 - 2 * jmax*vt) / (2 * jmax);
		val(1) = -(6 * a0 * a0 * amax - 3 * a0 * a0 * ap2 + 6 * amax * amax * ap2 - 6 * jmax * jmax * p0 + 6 * jmax * jmax * pt - 2 * a0 * a0 * a0 - 
			6 * amax * amax * amax - 6 * e12*jmax *jmax * v0 - 3 * amax*e12 *e12 * jmax * jmax + 6 * a0*jmax*v0 - 12 * amax*jmax*v0 + 6 * ap2*jmax*v0 +
			6 * ap2*jmax*vt + 3 * a0 * a0 * e12*jmax - 9 * amax * amax * e12*jmax + 6 * amax*ap2*e12*jmax) / (6 * jmax * jmax);
		return val;
	}

	MatrixX InverseJacobian(const VectorX & x) const
	{
		MatrixX j(2, 2);
		Real ap2 = x(0), e12 = x(1);
		Real det = (2 * a0 * a0 * ap2 - a0 * a0 * amax + 4 * amax*ap2*ap2 - 6 * amax*amax * ap2 + 2 * amax*amax*amax + 2 * amax*jmax*v0 -
			4 * ap2*jmax*v0 + 2 * amax*jmax*vt + 2 * amax*amax * e12*jmax - 4 * amax*ap2*e12*jmax) / (2 * jmax*jmax);
		j(0, 0) = ((2 * jmax*v0 - 2 * amax*ap2 - a0 * a0 + 3 * amax * amax + 2 * amax*e12*jmax) / (2 * jmax)) / det;
		j(1, 1) = (-(2 * ap2) / jmax) / det;
		j(0, 1) = -amax / det;
		j(1, 0) = (-3 * a0 * a0 + 6 * amax * amax + 6 * e12*jmax*amax + 6 * jmax*v0 + 6 * jmax*vt) / (6 * jmax * jmax) / det;
		return j;
	}
};

class case4 : public Function // PosTrapZeroNegTri ( a >= 0 ), output : ap2, e12, e34
{
public:
	Real a0, vt, v0, pt, p0, jmax, amax, vmax;

public:
	case4(const Real _a0, const Real _vt, const Real _v0, const Real _pt, const Real _p0, const Real _jmax, const Real _amax, const Real _vmax)
		: a0(_a0), vt(_vt), v0(_v0), pt(_pt), p0(_p0), jmax(_jmax), amax(_amax), vmax(_vmax) {}

	VectorX func(const VectorX& x) const
	{
		Real ap2 = x(0), e12 = x(1), e34 = x(2);
		VectorX val(3);
		val(0) = (-a0 * a0 + 2 * amax * amax + 2 * e12*jmax*amax - 2 * ap2 * ap2 + 2 * jmax*v0 - 2 * jmax*vt) / (2 * jmax);
		val(1) = (-a0 * a0 + 2 * amax * amax + 2 * e12*jmax*amax + 2 * jmax*v0 - 2 * jmax*vmax) / (2 * jmax);
		val(2) = (3 * a0 * a0 * ap2 - 6 * a0 * a0 * amax - 6 * amax * amax * ap2 + 6 * jmax * jmax * p0 - 6 * jmax * jmax * pt +
			2 * a0 * a0 * a0 + 6 * amax * amax * amax + 6 * e12*jmax * jmax * v0 + 6 * e34*jmax * jmax * v0 + 3 * amax*e12 * e12 * jmax * jmax -
			6 * a0*jmax*v0 + 12 * amax*jmax*v0 - 6 * ap2*jmax*v0 - 6 * ap2*jmax*vt - 3 * a0 * a0 * e12*jmax -
			3 * a0 * a0 * e34*jmax + 9 * amax * amax * e12*jmax + 6 * amax * amax * e34*jmax + 6 * amax*e12*e34*jmax * jmax - 6 * amax*ap2*e12*jmax) / (6 * jmax * jmax);
		return val;
	}

	MatrixX InverseJacobian(const VectorX & x) const
	{
		MatrixX j(3, 3);
		j.setZero();
		Real ap2 = x(0), e12 = x(1), e34 = x(2);
		j(0, 0) = -jmax / (2 * ap2);
		j(0, 1) = jmax / (2 * ap2);
		j(1, 1) = 1 / amax;
		j(2, 0) = -(-a0 * a0 + 2 * amax * amax + 2 * e12*jmax*amax + 2 * jmax*v0 + 2 * jmax*vt) / (2 * ap2*(-a0 * a0 + 2 * amax * amax + 2 * e12*jmax*amax + 2 * jmax*v0));
		j(2, 1) = (2 * a0 * a0 * ap2 - a0 * a0 * amax + 4 * amax*ap2 * ap2 - 6 * amax * amax * ap2 + 2 * amax * amax + 2 * amax*jmax*v0 - 
			4 * ap2*jmax*v0 + 2 * amax*jmax*vt + 2 * amax * amax * e12*jmax - 4 * amax*ap2*e12*jmax - 
			4 * amax*ap2*e34*jmax) / (2 * amax*ap2*(-a0 * a0 + 2 * amax * amax + 2 * e12*jmax*amax + 2 * jmax*v0));
		j(2, 2) = (2 * jmax) / (-a0 * a0 + 2 * amax * amax + 2 * e12*jmax*amax + 2 * jmax*v0);
		return j;
	}
};

class case5 : public Function // PosTrapNegTrap ( a >= 0 ), output : e12, e45
{
public:
	Real a0, vt, v0, pt, p0, jmax, amax, vmax;

public:
	case5(const Real _a0, const Real _vt, const Real _v0, const Real _pt, const Real _p0, const Real _jmax, const Real _amax, const Real _vmax)
		: a0(_a0), vt(_vt), v0(_v0), pt(_pt), p0(_p0), jmax(_jmax), amax(_amax), vmax(_vmax) {}

	VectorX func(const VectorX& x) const
	{
		Real e12 = x(0), e45 = x(1);
		VectorX val(2);

		val(0) = -(a0 * a0 - 2 * jmax*v0 + 2 * jmax*vt - 2 * amax*e12*jmax + 2 * amax*e45*jmax) / (2 * jmax);
		val(1) = (6 * jmax * jmax * p0 - 9 * a0 * a0 * amax - 6 * jmax * jmax * pt + 2 * a0 * a0 * a0 + 12 * amax * amax * amax + 
			6 * e12*jmax * jmax * v0 + 6 * e45*jmax * jmax * vt + 3 * amax*e12 * e12 * jmax * jmax + 3 * amax*e45 * e45 * jmax * jmax - 
			6 * a0*jmax*v0 + 18 * amax*jmax*v0 + 6 * amax*jmax*vt - 3 * a0 * a0 * e12*jmax + 15 * amax * amax * e12*jmax + 3 * amax * amax * e45*jmax) / (6 * jmax * jmax);

		return val;
	}

	MatrixX InverseJacobian(const VectorX & x) const
	{
		MatrixX j(2, 2);
		Real e12 = x(0), e45 = x(1);
		Real det = (amax*(2 * jmax*v0 + 2 * jmax*vt - a0 * a0 + 6 * amax * amax + 2 * amax*e12*jmax + 2 * amax*e45*jmax)) / (2 * jmax);
		j(0, 0) = ((amax * amax + 2 * e45*jmax*amax + 2 * jmax*vt) / (2 * jmax)) / det;
		j(1, 1) = (amax) / det;
		j(0, 1) = (amax) / det;
		j(1, 0) = -((-a0 * a0 + 5 * amax * amax + 2 * e12*jmax*amax + 2 * jmax*v0) / (2 * jmax)) / det;
		return j;
	}
};

Real makeRandLU(Real lower, Real upper);

int main()
{
	std::shared_ptr<Objective> obj = std::shared_ptr<Objective>(new Objective(6.25, 0, 10, 5, 15.4167, 0));
	NonlinearOptimization opt(2, 0, 0);
	opt.setObjectiveFunction(obj);
	VectorX initX(2), initX3(3);

	Real lower = -20, upper = 20;
	srand(time(NULL));
	Real value2 = rand();
	initX(0) = makeRandLU(0, 15);
	initX(1) = makeRandLU(-10, 0);
	cout << "--- initX ---" << endl;
	cout << "initX(0) : " << initX(0) << endl;
	cout << "initX(1) : " << initX(1) << endl;
	//initX(0) = 15; initX(1) = -6;
	initX3(0) = makeRandLU(0, 15);
	initX3(1) = makeRandLU(-10, 0);
	initX3(2) = makeRandLU(0, 2);
	cout << "--- initX3 ---" << endl;
	cout << "initX3(0) : " << initX3(0) << endl;
	cout << "initX3(1) : " << initX3(1) << endl;
	cout << "initX3(2) : " << initX3(2) << endl;

	cout << endl;
	clock_t time = clock();
	//for (int i = 0; i < 10000; i++)
		opt.solve(initX);
	cout << "computation time : " << clock() - time << endl;
	cout << opt.resultX << endl;
	cout << opt.resultFunc << endl;

	cout << "------- case 1 : PosTriNegTri -------" << endl;
	std::shared_ptr<case1> objfcn = std::shared_ptr<case1>(new case1(6.25, 0, 10, 5, 15.4167, 0));
	NewtonRaphson newton(2, 2);
	newton.setfunction(objfcn);
	time = clock();
	//for (int i = 0; i < 10000; i++)
		newton.solve(initX);
	cout << "computation time : " << clock() - time << endl;
	cout << newton.getResultX() << endl;

	cout << "------- case 2 : PosTriZeroNegTri -------" << endl;
	std::shared_ptr<case2> objfcn_case2 = std::shared_ptr<case2>(new case2(5, 11.25, 0, 94.7917, 0, 10, 0, 21.25));
	NewtonRaphson newtoncase2(3, 3);
	newtoncase2.setfunction(objfcn_case2);
	time = clock();
	newtoncase2.solve(initX3);
	cout << "computation time : " << clock() - time << endl;
	cout << newtoncase2.getResultX() << endl;
	//cout << objfcn_case2->func(newtoncase2.getResultX()) << endl;

	cout << "------- case 3 : PosTrapNegTri -------" << endl;
	initX(0) = makeRandLU(-10, 0);
	initX(1) = makeRandLU(0, 1);
	cout << "initX(0) : " << initX(0) << endl;
	cout << "initX(1) : " << initX(1) << endl;
	std::shared_ptr<case3> objfcn_case3 = std::shared_ptr<case3>(new case3(5, 11.25, 0, 28.5417, 0, 10, 10, 0));
	NewtonRaphson newtoncase3(2, 2);
	newtoncase3.setfunction(objfcn_case3);
	time = clock();
	newtoncase3.solve(initX);
	cout << "computation time : " << clock() - time << endl;
	cout << newtoncase3.getResultX() << endl;
	//cout << objfcn_case3->func(newtoncase3.getResultX()) << endl;

	cout << "------- case 4 : PosTrapZeroNegTri -------" << endl;
	initX3(0) = makeRandLU(-10, 0);
	initX3(1) = makeRandLU(0, 1);
	initX3(2) = makeRandLU(0, 1);
	//initX3(0) = -5; initX3(1) = 0.5; initX3(2) = 0.5;
	cout << "initX3(0) : " << initX3(0) << endl;
	cout << "initX3(1) : " << initX3(1) << endl;
	cout << "initX3(2) : " << initX3(2) << endl;
	std::shared_ptr<case4> objfcn_case4 = std::shared_ptr<case4>(new case4(5, 11.25, 0, 35.4167, 0, 10, 10, 13.75));
	NewtonRaphson newtoncase4(3, 3);
	newtoncase4.setfunction(objfcn_case4);
	time = clock();
	newtoncase4.solve(initX3);
	cout << "computation time : " << clock() - time << endl;
	cout << newtoncase4.getResultX() << endl;
	//cout << objfcn_case4->func(newtoncase4.getResultX()) << endl;

	cout << "------- case 5 : PosTrapNegTrap -------" << endl;
	initX(0) = makeRandLU(0, 2);
	initX(1) = makeRandLU(0, 2);
	//initX(0) = 0.5; initX(1) = 0.5;
	cout << "initX(0) : " << initX(0) << endl;
	cout << "initX(1) : " << initX(1) << endl;
	std::shared_ptr<case5> objfcn_case5 = std::shared_ptr<case5>(new case5(5, -1.25, 0, 31.6667, 0, 10, 10, 13.75));
	NewtonRaphson newtoncase5(2, 2);
	newtoncase5.setfunction(objfcn_case5);
	time = clock();
	newtoncase5.solve(initX);
	cout << "computation time : " << clock() - time << endl;
	cout << newtoncase5.getResultX() << endl;

	cout << "Program Complete" << endl;
	_getch();
	return 0;
}

Real makeRandLU(Real lower, Real upper)
{
	return lower + (upper - lower)*((double)rand() / 32767.0);
}