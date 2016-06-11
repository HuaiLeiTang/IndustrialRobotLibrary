#include <iostream>
#include <conio.h>
#include <rovin\Math\Function.h>
#include <memory>
#include <time.h>

#include "rovin\Optimizer\NonlinearOptimization.h"
#include "rovin\Optimizer\UnconstraintedOptimization.h"

using namespace std;
using namespace rovin;
using namespace irLib::Opt;

Real makeRandLU(Real lower, Real upper);

class Step2PosTriZeroNegTri : public Function // unknown : ap1, ap2
{
public:
	Real tsyn, a0, vt, v0, pt, p0, jmax, amax, vmax;

public:
	Step2PosTriZeroNegTri(const Real _tsyn, const Real _a0, const Real _vt, const Real _v0, const Real _pt, const Real _p0, const Real _jmax, const Real _amax, const Real _vmax)
		: tsyn(_tsyn), a0(_a0), vt(_vt), v0(_v0), pt(_pt), p0(_p0), jmax(_jmax), amax(_amax), vmax(_vmax) {}

	VectorX func(const VectorX& x) const
	{
		//VectorX val(3);
		//Real ap1 = x(0), ap2 = x(1), t34 = x(2);
		//val(0) = -(a0 - 2 * ap1 + 2 * ap2 - jmax*t34 + jmax*tsyn) / jmax;
		//val(1) = -(a0 * a0 - 2 * ap1 * ap1 + 2 * ap2 * ap2 - 2 * jmax*v0 + 2 * jmax*vt) / (2 * jmax);
		//val(2) = (6 * a0 * a0 * ap2 - 6 * a0 * a0 * ap1 - 12 * ap1 * ap1 * ap2 + 6 * jmax * jmax * p0 - 6 * jmax * jmax * pt + 
		//	2 * a0 * a0 * a0 + 6 * ap1 * ap1 * ap1 + 6 * ap2 * ap2 * ap2 + 6 * jmax * jmax * t34*v0 - 6 * a0*jmax*v0 + 12 * ap1*jmax*v0 - 
		//	12 * ap2*jmax*v0 - 3 * a0 * a0 * jmax*t34 + 6 * ap1 * ap1 * jmax*t34) / (6 * jmax * jmax);

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

class Step2PosTrapZeroNegTri : public Function // unknown : ap2, t23
{
public:
	Real tsyn, a0, vt, v0, pt, p0, jmax, amax, vmax;

public:
	Step2PosTrapZeroNegTri(const Real _tsyn, const Real _a0, const Real _vt, const Real _v0, const Real _pt, const Real _p0, const Real _jmax, const Real _amax, const Real _vmax)
		: tsyn(_tsyn), a0(_a0), vt(_vt), v0(_v0), pt(_pt), p0(_p0), jmax(_jmax), amax(_amax), vmax(_vmax) {}

	VectorX func(const VectorX& x) const
	{
		VectorX val(2);
		Real ap2 = x(0), t23 = x(1);
		val(0) = (-a0*a0 + 2 * amax*amax + 2 * jmax*t23*amax - 2 * ap2*ap2 + 2 * jmax*v0 - 2 * jmax*vt) / (2 * jmax);
		val(1) = (6 * a0*amax*amax + 6 * jmax*jmax * p0 - 6 * jmax*jmax * pt - a0*a0*a0 - 6 * amax*amax*amax + 6 * ap2*ap2*ap2 + 6 * jmax*jmax * tsyn*v0 - 
			3 * amax*jmax*jmax * t23*t23 - 9 * amax*amax * jmax*t23 - 3 * a0*a0 * jmax*tsyn + 6 * amax*amax * jmax*tsyn + 6 * a0*amax*jmax*t23 + 
			6 * amax*jmax*jmax * t23*tsyn) / (6 * jmax * jmax);
		return val;
	}

	MatrixX InverseJacobian(const VectorX & x) const
	{
		MatrixX j(2, 2);
		Real ap2 = x(0), t23 = x(1);
		Real det = -(amax*ap2*(2 * a0 - 3 * amax + 3 * ap2 - 2 * jmax*t23 + 2 * jmax*tsyn)) / (jmax*jmax);
		j(0, 0) = ((amax*(2 * a0 - 3 * amax - 2 * jmax*t23 + 2 * jmax*tsyn)) / (2 * jmax)) / det;
		j(1, 1) = (-(2 * ap2) / jmax) / det;
		j(0, 1) = -(amax) / det;
		j(1, 0) = -((3 * ap2 * ap2) / (jmax * jmax)) / det;
		return j;
	}
};

class Step2PosTrapZeroNegTrap : public Function // unknown : t23, t67
{
public:
	Real tsyn, a0, vt, v0, pt, p0, jmax, amax, vmax;

public:
	Step2PosTrapZeroNegTrap(const Real _tsyn, const Real _a0, const Real _vt, const Real _v0, const Real _pt, const Real _p0, const Real _jmax, const Real _amax, const Real _vmax)
		: tsyn(_tsyn), a0(_a0), vt(_vt), v0(_v0), pt(_pt), p0(_p0), jmax(_jmax), amax(_amax), vmax(_vmax) {}

	VectorX func(const VectorX& x) const
	{
		VectorX val(2);
		Real t23 = x(0), t67 = x(1);
		val(0) = -(a0 * a0 - 2 * jmax*v0 + 2 * jmax*vt - 2 * amax*jmax*t23 + 2 * amax*jmax*t67) / (2 * jmax);
		val(1) = -(6 * jmax * jmax * pt - 6 * jmax * jmax * p0 - 6 * a0*amax *amax + a0 * a0 * a0 + 12 * amax * amax*amax - 6 * jmax * jmax * tsyn*v0 +
			3 * amax*jmax * jmax * t23 * t23 + 3 * amax*jmax * jmax * t67 * t67 + 9 * amax * amax * jmax*t23 + 9 * amax * amax * jmax*t67 +
			3 * a0 * a0 * jmax*tsyn - 6 * amax * amax * jmax*tsyn - 6 * a0*amax*jmax*t23 - 6 * amax*jmax * jmax * t23*tsyn) / (6 * jmax * jmax);
		return val;
	}

	MatrixX InverseJacobian(const VectorX & x) const
	{
		MatrixX j(2, 2);
		Real t23 = x(0), t67 = x(1);
		Real det = -(amax * amax * (3 * amax - a0 + jmax*t23 + jmax*t67 - jmax*tsyn)) / jmax;
		j(0, 0) = (-(amax*(3 * amax + 2 * jmax*t67)) / (2 * jmax)) / det;
		j(1, 1) = (amax) / det;
		j(0, 1) = -(-amax) / det;
		j(1, 0) = -((amax*(2 * a0 - 3 * amax - 2 * jmax*t23 + 2 * jmax*tsyn)) / (2 * jmax)) / det;
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
	clock_t time;

	// PosTriZeroNegTri Case Experiment
	cout << "------- case 1 : PosTriZeroNegTri -------" << endl;
	std::shared_ptr<Step2PosTriZeroNegTri> obj1 = std::shared_ptr<Step2PosTriZeroNegTri>(new Step2PosTriZeroNegTri(3, 5, 6.25, 0, 19.7917, 0, 10, 100, 100));
	vec2(0) = makeRandLU(0, 15); vec2(1) = makeRandLU(-10, 0);
	cout << "initial ap1 : " << vec2(0) << endl;
	cout << "initial ap2 : " << vec2(1) << endl;
	solver.setxN(2); solver.setfN(2); solver.setfunction(obj1);
	time = clock();
	solver.solve(vec2);
	cout << "Computation time : " << clock() - time << endl;
	cout << "Result" << endl << solver.getResultX() << endl << endl;

	// PosTrapZeroNegTri Case Experiment
	cout << "------- case 2 : PosTrapZeroNegTri -------" << endl;
	std::shared_ptr<Step2PosTrapZeroNegTri> obj2 = std::shared_ptr<Step2PosTrapZeroNegTri>(new Step2PosTrapZeroNegTri(6, 5, 18.75, 0, 124.7917, 0, 10, 15, 100));
	vec2(0) = makeRandLU(-15, 0); vec2(1) = makeRandLU(0, 1);
	cout << "initial ap2 : " << vec2(0) << endl;
	cout << "initial t23 : " << vec2(1) << endl;
	solver.setxN(2); solver.setfN(2); solver.setfunction(obj2);
	time = clock();
	solver.solve(vec2);
	cout << "Computation time : " << clock() - time << endl;
	cout << "Result" << endl << solver.getResultX() << endl << endl;

	// PosTrapZeroNegTrap Case Experiment
	cout << "------- case 3 : PosTrapZeroNegTrap -------" << endl;
	std::shared_ptr<Step2PosTrapZeroNegTrap> obj3 = std::shared_ptr<Step2PosTrapZeroNegTrap>(new Step2PosTrapZeroNegTrap(6, 5, -6.25, 0, 41.0417, 0, 10, 10, 100));
	vec2(0) = makeRandLU(0, 1); vec2(1) = makeRandLU(0, 1);
	cout << "initial t23 : " << vec2(0) << endl;
	cout << "initial t67 : " << vec2(1) << endl;
	solver.setxN(2); solver.setfN(2); solver.setfunction(obj3);
	time = clock();
	solver.solve(vec2);
	cout << "Computation time : " << clock() - time << endl;
	cout << "Result" << endl << solver.getResultX() << endl << endl;

	cout << "Program Complete" << endl;
	_getch();
	return 0;
}

Real makeRandLU(Real lower, Real upper)
{
	return lower + (upper - lower)*((double)rand() / 32767.0);
}