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

typedef struct datastr
{
	Real vt;
	Real v0;
	Real jmax;
	Real a0;
	Real pt;
	Real p0;
} datastr;

void case2fcn(const VectorX& x, VectorX& f, MatrixX& ig, void* f_data)
{
	Real vt, v0, jmax, a0, pt, p0;
	datastr* d = reinterpret_cast<datastr*>(f_data);
	vt = d->vt;
	v0 = d->v0;
	jmax = d->jmax;
	a0 = d->a0;
	pt = d->pt;
	p0 = d->p0;

	Real ap1 = x(0), ap2 = x(1);
	f(0) = vt - v0 - (ap1*ap1 / jmax - 0.5*a0*a0 / jmax - ap2*ap2 / jmax);
	f(1) = pt - p0 - (-2 * vt*ap2 / jmax - ap2*ap2*ap2 / (jmax*jmax) + 2 * v0*ap1 / jmax - v0*a0 / jmax + a0*a0*a0 / (3 * jmax*jmax) +
		ap1*ap1*ap1 / (jmax*jmax) - a0*a0*ap1 / (jmax*jmax));

	Real det = -(2 * (a0*a0 * ap2 - 3 * ap1*ap1 * ap2 + 3 * ap1*ap2*ap2 + 2 * jmax*vt*ap1 - 2 * jmax*v0*ap2)) / (jmax*jmax*jmax);
	ig(0, 0) = ((2 * vt) / jmax + (3 * ap2 * ap2) / (jmax*jmax)) / det;
	ig(0, 1) = -((2 * ap2) / jmax) / det;
	ig(1, 0) = -(a0 * a0 / (jmax *jmax) - (2 * v0) / jmax - (3 * ap1 * ap1) / (jmax*jmax)) / det;
	ig(1, 1) = (-(2 * ap1) / jmax) / det;
}

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
		val(1) = pt - p0 - (-2 * vt*ap2 / jmax - ap2*ap2*ap2 / (jmax*jmax) + 2 * v0*ap1 / jmax - v0*a0 / jmax + a0*a0*a0 / (3 * jmax*jmax) + 
			ap1*ap1*ap1 / (jmax*jmax) - a0*a0*ap1 / (jmax*jmax));
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

class case2 : public Function // PosTriZeroNegTri ( a >= 0 ), output : ap1, ap2, e34 , there exist closed-form
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

	void closedform()
	{
		Real ap1, ap2, e23;
		ap1 = 1 / std::sqrt(2) * std::sqrt(a0*a0 + 2*jmax*(vmax - v0));
		ap2 = -1 / std::sqrt(2) * std::sqrt(-a0*a0 + 2*ap1*ap1 + 2*jmax*v0 - 2*jmax*vt);
		e23 = (6 * a0 * a0 * ap1 - 3 * a0 * a0 * ap2 + 6 * ap1 * ap1 * ap2 - 6 * jmax * jmax * p0 + 6 * jmax * jmax * pt -
			2 * a0 * a0 * a0 - 6 * ap1 * ap1 * ap1 + 6 * a0*jmax*v0 - 12 * ap1*jmax*v0 +
			6 * ap2 * jmax*v0 + 6 * ap2 * jmax*vt) / (6 * jmax * jmax * v0 - 3 * a0 * a0 * jmax + 6 * ap1 * ap1 * jmax);
			 
		cout << "ap1 : " << ap1 << endl;
		cout << "ap2 : " << ap2 << endl;
		cout << "e23 : " << e23 << endl;
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

	void closedform()
	{
		Real ap2, e12;
		Real a, b, c, d, e;
		a = 12;
		b = -24 * amax;
		c = 12 * amax * amax + 24 * jmax*vt;
		d = -48 * amax*jmax*vt;
		e = 8 * std::pow(a0, 3) * amax - 3 * std::pow(a0, 4) - 6 * a0 * a0 * amax * amax - 12 * jmax * jmax * v0 * v0 + 12 * jmax * jmax * vt * vt + 
			12 * a0 * a0 * jmax*v0 + 12 * amax * amax * jmax*v0 + 12 * amax * amax * jmax*vt + 24 * amax*jmax * jmax * p0 - 
			24 * amax*jmax * jmax * pt - 24 * a0*amax*jmax*v0;

		Real p, q, det0, det1, Q, S;
		p = (8 * a*c - 3 * b*b) / (8 * a*a);
		q = (b*b*b - 4 * a*b*c + 8 * a*a*d) / (8 * a*a*a);
		det0 = c*c - 3 * b*d + 12 * a*e;
		det1 = 2 * c*c*c - 9 * b*c*d + 27 * b*b*e + 27 * a*d*d - 72 * a*c*e;
		Q = std::pow(((det1 + std::sqrt(det1*det1 - 4 * det0*det0*det0)) / 2), (1.0 / 3.0));
		S = 1.0 / 2.0* std::sqrt(-2.0 / 3.0*p + 1.0 / 3.0 / a*(Q + det0 / Q));

		ap2 = -b / 4 / a + S - 0.5*sqrt(-4 * S*S - 2 * p - q / S);
		e12 = (a0 * a0 - 2 * amax * amax + 2 * ap2 * ap2 - 2 * jmax*v0 + 2 * jmax*vt) / (2 * amax*jmax);

		cout << "ap2 : " << ap2 << endl;
		cout << "e12 : " << e12 << endl;
	}
};

class case4 : public Function // PosTrapZeroNegTri ( a >= 0 ), output : ap2, e12, e34, there exist closed-form.
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

	void closedform()
	{
		Real ap2, e12, e34;
		
		e12 = (-a0 * a0 + 2 * amax * amax + 2 * jmax*v0 - 2 * jmax*vmax) / (-2*jmax*amax);
		ap2 = -1 / std::sqrt(2) * std::sqrt(-a0 * a0 + 2 * amax * amax + 2 * e12*jmax*amax + 2 * jmax*v0 - 2 * jmax*vt);
		e34 = ( 3 * a0 * a0 * ap2 - 6 * a0 * a0 * amax - 6 * amax * amax * ap2 + 6 * jmax * jmax * p0 - 6 * jmax * jmax * pt +
			2 * a0 * a0 * a0 + 6 * amax * amax * amax + 6 * e12*jmax * jmax * v0 + 3 * amax*e12 * e12 * jmax * jmax -
			6 * a0*jmax*v0 + 12 * amax*jmax*v0 - 6 * ap2*jmax*v0 - 6 * ap2*jmax*vt - 3 * a0 * a0 * e12*jmax
			+ 9 * amax * amax * e12*jmax - 6 * amax*ap2*e12*jmax ) / (-6 * jmax * jmax * v0 + 3 * a0 * a0 *jmax - 6 * amax * amax *jmax - 6 * amax*e12*jmax* jmax);
		 
		cout << "ap2 : " << ap2 << endl;
		cout << "e12 : " << e12 << endl;
		cout << "e34 : " << e34 << endl;
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

	void closedform()
	{
		Real e12, e45;
		e45 = -(2 * jmax*vt + 3 * amax * amax - std::sqrt((std::pow(a0, 4)) / 2 - (4 * std::pow(a0, 3) * amax) / 3 + 
			std::pow(amax, 4) + a0*a0 * amax*amax + 2 * jmax*jmax * v0*v0 + 2 * jmax*jmax * vt*vt - 2 * a0*a0 * jmax*v0 - 
			2 * amax*amax * jmax*v0 - 2 * amax * amax * jmax*vt - 4 * amax*jmax*jmax * p0 + 4 * amax*jmax * jmax * pt + 4 * a0*amax*jmax*v0)) / (2 * amax*jmax);
		e12 = (a0 * a0 - 2 * jmax*v0 + 2 * jmax*vt + 2 * amax*e45*jmax) / (2 * amax*jmax);
		cout << "e12 : " << e12 << endl;
		cout << "e34 : " << e45 << endl;
	}
};

class case6
{
public:
	Real a0, vt, v0, pt, p0, jmax, amax, vmax;

public:
	case6(const Real _a0, const Real _vt, const Real _v0, const Real _pt, const Real _p0, const Real _jmax, const Real _amax, const Real _vmax)
		: a0(_a0), vt(_vt), v0(_v0), pt(_pt), p0(_p0), jmax(_jmax), amax(_amax), vmax(_vmax) {}
	
	void closedform()
	{
		Real e12, e34, e56;
		e12 = -(-a0 * a0 + 2 * amax * amax + 2 * jmax*v0 - 2 * jmax*vmax) / (2 * amax*jmax);
		e56 = -(amax * amax - jmax*vmax + jmax*vt) / (amax*jmax);
		e34 = -(8 * std::pow(a0, 3) * amax - 3 * std::pow(a0, 4) - 6 * a0*a0 * amax*amax - 12 * jmax*jmax * v0*v0 + 24 * jmax*jmax * vmax*vmax -
			12 * jmax*jmax * vt*vt + 12 * a0*a0 * jmax*v0 + 12 * amax*amax * jmax*v0 + 24 * amax*amax * jmax*vmax +
			12 * amax*amax * jmax*vt + 24 * amax*jmax*jmax * p0 - 24 * amax*jmax*jmax * pt - 24 * a0*amax*jmax*v0) / (24 * amax*jmax*jmax * vmax);
		
		cout << "e12 : " << e12 << endl;
		cout << "e34 : " << e34 << endl;
		cout << "e56 : " << e56 << endl;
	}


};

class case7
{
public:
	Real a0, vt, v0, pt, p0, jmax, amax, vmax;
	case7(const Real _a0, const Real _vt, const Real _v0, const Real _pt, const Real _p0, const Real _jmax, const Real _amax, const Real _vmax)
		: a0(_a0), vt(_vt), v0(_v0), pt(_pt), p0(_p0), jmax(_jmax), amax(_amax), vmax(_vmax) {}

	void closedform()
	{
		Real ap1, e34;

		Real a, b, c, d, e;
		a = 12;
		b = 24 * amax;
		c = -12 * a0*a0 + 12 * amax*amax + 24 * jmax*v0;
		d = -24 * a0*a0 * amax + 48 * amax*jmax*v0;
		e = 8 * std::pow(a0, 3) * amax + 3 * std::pow(a0, 4) - 6 * a0*a0 * amax*amax + 12 * jmax*jmax * v0*v0 - 12 * jmax*jmax * vt*vt
			- 12 * a0 * a0 * jmax*v0 + 12 * amax * amax * jmax*v0 + 12 * amax * amax * jmax*vt + 24 * amax*jmax * jmax * p0
			- 24 * amax*jmax*jmax * pt - 24 * a0*amax*jmax*v0;

		Real p, q, det0, det1, Q, S;
		p = (8 * a*c - 3 * b*b) / (8 * a*a);
		q = (b*b*b - 4 * a*b*c + 8 * a*a*d) / (8 * a*a*a);
		det0 = c*c - 3 * b*d + 12 * a*e;
		det1 = 2 * c*c*c - 9 * b*c*d + 27 * b*b*e + 27 * a*d*d - 72 * a*c*e;
		Q = std::pow(((det1 + std::sqrt(det1*det1 - 4 * det0*det0*det0)) / 2), (1.0 / 3.0));
		S = 1.0 / 2.0* std::sqrt(-2.0 / 3.0*p + 1.0 / 3.0 / a*(Q + det0 / Q));

		ap1 = -b / 4 / a + S + 0.5*sqrt(-4 * S*S - 2 * p - q / S);
		e34 = -(a0 * a0 + 2 * amax * amax - 2 * ap1 * ap1 - 2 * jmax*v0 + 2 * jmax*vt) / (2 * amax*jmax);

		cout << "ap1 : " << ap1 << endl;
		cout << "e34 : " << e34 << endl;
	}

};

class case8
{
public:
	Real a0, vt, v0, pt, p0, jmax, amax, vmax;
	case8(const Real _a0, const Real _vt, const Real _v0, const Real _pt, const Real _p0, const Real _jmax, const Real _amax, const Real _vmax)
		: a0(_a0), vt(_vt), v0(_v0), pt(_pt), p0(_p0), jmax(_jmax), amax(_amax), vmax(_vmax) {}

	void closedform()
	{
		Real ap1, e23, e45;
		ap1 = std::sqrt(0.5*a0*a0 + vmax*jmax - v0*jmax);

		Real e01, e12, e34, e56;
		e01 = (ap1 - a0) / jmax;
		e12 = ap1 / jmax;
		e34 = amax / jmax;
		e56 = e34;

		Real v01, v12, v34, v45, v56;
		v01 = 0.5*(a0 + ap1)*e01;
		v12 = 0.5*ap1*e12;
		v34 = -0.5*amax*e34;
		v56 = v34;

		e45 = (v0 - vt + v01 + v12 + v34 + v56) / amax;

		v45 = -amax*e45;

		Real xi01, xi12, xi23, xi34, xi45, xi56;
		xi01 = v0*e01 + 0.5*a0*e01*e01 + 1.0 / 6.0*jmax*e01*e01*e01;
		xi12 = (v0 + v01)*e12 + 0.5*ap1*e12*e12 - 1.0 / 6.0*jmax*e12*e12*e12;
		xi34 = (v0 + v01 + v12)*e34 - 1.0 / 6.0 *jmax*e34*e34*e34;
		xi45 = (vt - v45 - v56)*e45 - 0.5*amax*e45*e45;
		xi56 = (vt - v56)*e56 - 0.5*amax*e56*e56 + 1.0 / 6.0*jmax*e56*e56*e56;

		e23 = (pt - p0 - xi01 - xi12 - xi34 - xi45 - xi56) / vmax;

		cout << "ap1 : " << ap1 << endl;
		cout << "e23 : " << e23 << endl;
		cout << "e45 : " << e45 << endl;
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

	////////////////////////
	cout << "------- NewtonRaphsonfcn Test -------" << endl;
	datastr* d = new datastr;
	d->a0 = 5;
	d->jmax = 10;
	d->p0 = 0;
	d->pt = 15.4167;
	d->v0 = 0;
	d->vt = 6.25;

	NewtonRaphsonfcn testNewton(2, 2);
	testNewton.settolerance(1E-10);
	testNewton.setfunction(case2fcn);
	testNewton.setfdata(d);
	time = clock();
	testNewton.solve(initX);
	cout << "computation time : " << clock() - time << endl;
	cout << testNewton.getResultX() << endl;
	delete d;
	////////////////////////

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
	cout << "closedform" << endl;
	objfcn_case2->closedform();
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
	objfcn_case3->closedform();

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
	cout << "closedform" << endl;
	objfcn_case4->closedform();
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
	cout << "closedform" << endl;
	objfcn_case5->closedform();

	cout << "------- case 6 : PosTrapZeroNegTrap -------" << endl;
	case6 objfcn_case6(5, -6.25, 0, 41.0417, 0, 10, 10, 13.75);
	objfcn_case6.closedform();


	cout << "------- case 7 : PosTriNegTrap -------" << endl;
	case7 objfcn_case7(5, -28.75, 0, -32.0833, 0, 10, 15, 13.75);
	objfcn_case7.closedform();

	cout << "------- case 8 : PosTriZeroNegTrap -------" << endl;
	case8 objcn_case8(5, -28.75, 0, -27.7083, 0, 10, 15, 8.75);
	objcn_case8.closedform();

	cout << "Program Complete" << endl;
	_getch();
	return 0;
}

Real makeRandLU(Real lower, Real upper)
{
	return lower + (upper - lower)*((double)rand() / 32767.0);
}