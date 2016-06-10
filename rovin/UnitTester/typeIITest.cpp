#include <iostream>
#include <conio.h>
#include <time.h>
#include "rovin\Math\Function.h"
#include "rovin\Optimizer\NonlinearOptimization.h"

using namespace std;
using namespace rovin;

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
		f1 = vt - v0 - (ap1*ap1 / jmax + 0.5*(a0 - 1)*ap1 / jmax - 0.5*a0*a0 / jmax - ap2*ap2 / jmax);
		f2 = pt - p0 - (-2 * vt*ap2 / jmax - 2 * ap2*ap2*ap2 / (3 * jmax*jmax) + ap1*ap1*ap1 / (jmax*jmax) + 2 * v0*ap1 / jmax - a0*a0*ap1 / (jmax*jmax) - v0*a0 / jmax + a0*a0*a0 / (3 * jmax*jmax));
		val(0) = f1*f1 + f2*f2;
		return val;
	}

	void result(Real ap1, Real ap2)
	{
		Real f1 = 0, f2 = 0;
		f1 = ap1*ap1 / jmax + 0.5*(a0 - 1)*ap1 / jmax - 0.5*a0*a0 / jmax - ap2*ap2 / jmax;
		f2 = -2 * vt*ap2 / jmax - 2 * ap2*ap2*ap2 / (3 * jmax*jmax) + ap1*ap1*ap1 / (jmax*jmax) + 2 * v0*ap1 / jmax - a0*a0*ap1 / (jmax*jmax) - v0*a0 / jmax + a0*a0*a0 / (3 * jmax*jmax);
		cout << "vt - v0 : " << vt - v0 << endl;
		cout << "pt - p0 : " << pt - p0 << endl;
		cout << "f1 : " << f1 << endl;
		cout << "f2 : " << f2 << endl;
	}
};

int main()
{
	std::shared_ptr<Objective> obj = std::shared_ptr<Objective>(new Objective(6.25, 0, 10, 5, 9.1667, 0));
	NonlinearOptimization opt(2, 0, 0);
	opt.setObjectiveFunction(obj);
	VectorX initX(2);
	initX(0) = 9;
	initX(1) = -4;
	clock_t time = clock();
	for (int i = 0; i < 1000; i++)
		opt.solve(initX);
	cout << "computation time : " << clock() - time << endl;

	cout << opt.resultX << endl;
	cout << opt.resultFunc << endl;



	cout << "Program Complete" << endl;
	_getch();
	return 0;
}