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

class ComplexNumber
{
public:
	Real _rN; // real number
	Real _iN; // imagine number

public:
	ComplexNumber(Real rN = 0, Real iN = 0) : _rN(rN), _iN(iN) {}

	friend ComplexNumber operator+(const ComplexNumber& cn, const Real num)
	{
		ComplexNumber tmpCN;
		tmpCN._rN = cn._rN + num;
		tmpCN._iN = cn._iN;
		return tmpCN;
	}

	friend ComplexNumber operator+(const Real num, const ComplexNumber& cn)
	{
		ComplexNumber tmpCN;
		tmpCN._rN = cn._rN + num;
		tmpCN._iN = cn._iN;
		return tmpCN;
	}

	friend ComplexNumber operator+(const ComplexNumber& cn1, const ComplexNumber& cn2)
	{
		ComplexNumber tmpCN;
		tmpCN._rN = cn1._rN + cn2._rN;
		tmpCN._iN = cn1._iN + cn2._iN;
		return tmpCN;
	}

	friend ComplexNumber operator-(const ComplexNumber& cn, const Real num)
	{
		ComplexNumber tmpCN;
		tmpCN._rN = cn._rN - num;
		tmpCN._iN = cn._iN;
		return tmpCN;
	}

	friend ComplexNumber operator-(const Real num, const ComplexNumber& cn)
	{
		ComplexNumber tmpCN;
		tmpCN._rN = num - cn._rN;
		tmpCN._iN = cn._iN;
		return tmpCN;
	}

	friend ComplexNumber operator-(const ComplexNumber& cn1, const ComplexNumber& cn2)
	{
		ComplexNumber tmpCN;
		tmpCN._rN = cn1._rN - cn2._rN;
		tmpCN._iN = cn1._iN - cn2._iN;
		return tmpCN;
	}

	friend ComplexNumber operator*(const ComplexNumber& cn, const Real num)
	{
		ComplexNumber tmpCN;
		tmpCN._rN = cn._rN * num;
		tmpCN._iN = cn._iN * num;
		return tmpCN;
	}

	friend ComplexNumber operator*(const Real num, const ComplexNumber& cn)
	{
		ComplexNumber tmpCN;
		tmpCN._rN = cn._rN * num;
		tmpCN._iN = cn._iN * num;
		return tmpCN;
	}

	friend ComplexNumber operator*(const ComplexNumber& cn1, const ComplexNumber& cn2)
	{
		ComplexNumber tmpCN;
		Real a, b, c, d;
		a = cn1._rN; b = cn1._iN;
		c = cn2._rN; d = cn2._iN;
		tmpCN._rN = a*c - b*d;
		tmpCN._iN = a*d + b*c;
		return tmpCN;
	}

	friend ComplexNumber operator/(const Real num, const ComplexNumber& cn)
	{
		Real a, b;
		a = cn._rN; b = cn._iN;
		ComplexNumber tmpCN(a, b);
		tmpCN._rN = num*a / (a*a + b*b);
		tmpCN._iN = -num*b / (a*a + b*b);
		return tmpCN;
	}

	friend ComplexNumber operator/(const ComplexNumber& cn, const Real num)
	{
		ComplexNumber tmpCN(cn._rN, cn._iN);
		tmpCN._rN /= num;
		tmpCN._iN /= num;
		return tmpCN;
	}

	friend ostream& operator<<(ostream& os, const ComplexNumber& cn)
	{
		if (cn._iN == 0)
		{
			os << cn._rN << "+0.0i";
		}
		else if (cn._iN > 0)
		{
			os << cn._rN << " + " << cn._iN << "i";
		}
		else
		{
			os << cn._rN << " - " << std::abs(cn._iN) << "i";
		}
		return os;
	}

	void sqrt()
	{
		Real r, w, th;
		r = std::sqrt(_rN*_rN + _iN*_iN);
		w = std::sqrt(r);
		th = atan(_iN / _rN);
		_rN = w*cos(th / 2);
		_iN = w*sin(th / 2);
	}

	static ComplexNumber sqrt(const ComplexNumber& cn)
	{
		ComplexNumber tmpCN;
		Real r, w, th;
		Real ep = 1E-8;
		if (abs(cn._rN) > ep && abs(cn._iN) > ep)
		{
			r = std::sqrt(cn._rN*cn._rN + cn._iN*cn._iN);
			w = std::sqrt(r);
			th = atan(cn._iN / cn._rN);
			tmpCN._rN = w*cos(th / 2);
			tmpCN._iN = w*sin(th / 2);
		}
		else if (abs(cn._iN) < ep) // pure real
		{
			if (cn._rN >= 0)
				tmpCN._rN = std::sqrt(cn._rN);
			else
				tmpCN._iN = std::sqrt(abs(cn._rN));
		}
		else if (abs(cn._rN) < ep) // pure imaginary
		{
			tmpCN._rN = std::sqrt(std::abs(cn._iN)) / std::sqrt(2);
			tmpCN._iN = tmpCN._rN;
			if (cn._iN < 0)
				tmpCN._rN = -tmpCN._rN;
		}
		else // zero
		{
			tmpCN._rN = 0;
			tmpCN._iN = 0;
		}

		return tmpCN;
	}

	void cuberoot()
	{
		int n = 3;
		Real r, w, th;
		Real ep = 1E-8;
		if (abs(_rN) > ep && abs(_iN) > ep)
		{
			r = std::sqrt(_rN*_rN + _iN*_iN);
			w = std::pow(r, 1.0 / (double)(n));
			th = atan(_iN / _rN);
			_rN = w*cos(th / n);
			_iN = w*sin(th / n);
		}
		else if (abs(_iN) < ep) // pure real
		{
			if (_rN >= 0)
				_rN = std::pow(_rN, 1.0 / 3.0);
			else
			{
				_rN = std::pow(std::abs(_rN), 1.0 / 3.0) * 0.5;
				_iN = std::pow(std::abs(_rN), 1.0 / 3.0) * 0.866025403784439;
			}
		}
		else if (abs(_rN) < ep) // pure imaginary
		{
			if (_iN > 0)
			{
				_rN = std::pow(_rN, 1.0 / 3.0) * 0.866025403784439;
				_iN = std::pow(_rN, 1.0 / 3.0) * 0.5;
			}
			else
				_iN = std::pow(std::abs(_iN), 1.0 / 3.0);
		}
		else // zero
		{
			_rN = 0;
			_iN = 0;
		}
	}

	void root(const int n)
	{
		Real r, w, th;
		r = std::sqrt(_rN*_rN + _iN*_iN);
		w = std::pow(r, 1.0 / (double)(n));
		th = atan(_iN / _rN);
		_rN = w*cos(th / n);
		_iN = w*sin(th / n);
	}
};

void calcQuarticAlgebraicEqn(Real a, Real b, Real c, Real d, Real e, ComplexNumber* x)
{
	Real p, q, det0, det1;

	p = (8 * a*c - 3 * b*b) / (8 * a*a);
	q = (b*b*b - 4 * a*b*c + 8 * a*a*d) / (8 * a*a*a);
	det0 = c*c - 3 * b*d + 12 * a*e;
	det1 = 2 * c*c*c - 9 * b*c*d + 27 * b*b*e + 27 * a*d*d - 72 * a*c*e;

	Real det = det1*det1 - 4 * det0*det0*det0;
	ComplexNumber Q, S;

	if (det >= 0)
	{
		Q._rN = (det1 + std::sqrt(std::abs(det))) / 2;
	}
	else
	{
		Q._rN = det1 / 2;
		Q._iN = std::sqrt(std::abs(det)) / 2;
	}

	//Q.root(3);
	Q.cuberoot();

	//cout << "Q : " << Q << endl;

	S = 1.0 / 2.0 * ComplexNumber::sqrt(-2.0 / 3.0*p + 1.0 / 3.0 / a*(Q + det0 / Q));

	//cout << "S : " << S << endl;
	//cout << "-4 * S*S - 2 * p + q / S : " << -4 * S*S - 2 * p + q / S << endl;
	//cout << "-4 * S*S - 2 * p - q / S : " << -4 * S*S - 2 * p - q / S << endl;

	x[0] = -(b / 4 / a) - S + 0.5*ComplexNumber::sqrt(-4 * S*S - 2 * p + q / S);
	x[1] = -(b / 4 / a) - S - 0.5*ComplexNumber::sqrt(-4 * S*S - 2 * p + q / S);
	x[2] = -(b / 4 / a) + S + 0.5*ComplexNumber::sqrt(-4 * S*S - 2 * p - q / S);
	x[3] = -(b / 4 / a) + S - 0.5*ComplexNumber::sqrt(-4 * S*S - 2 * p - q / S);
}

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

//////////////////////// case 2 ////////////////////////
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

	void closedform()
	{
		Real a, b, c, d, e;
		a = 12;
		b = -24 * amax;
		c = 12 * a0*a0 + 12 * amax * amax + 24 * jmax*vt - 24 * amax*jmax*tsyn - 
			24 * a0*amax - 24 * jmax*v0;
		d = 0;
		e = 3 * std::pow(a0, 4) - 8 * std::pow(a0, 3) * amax + 6 * a0*a0 * amax * amax + 12 * jmax*jmax * v0*v0 + 12 * jmax*jmax * vt*vt - 
			12 * a0*a0 * jmax*v0 - 12 * amax * amax * jmax*v0 + 12 * a0*a0 * jmax*vt + 12 * amax * amax * jmax*vt - 
			24 * jmax*jmax * v0*vt - 24 * amax*jmax*jmax * p0 + 24 * amax*jmax*jmax * pt + 24 * a0*amax*jmax*v0 - 
			24 * a0*amax*jmax*vt - 24 * amax*jmax*jmax * tsyn*vt;

		ComplexNumber x[4];
		calcQuarticAlgebraicEqn(a, b, c, d, e, x);
		
		Real ep = 1E-8, sol = -RealMax;
		for (int i = 0; i < 4; i++)
		{
			if (abs(x[i]._iN) < ep)
			{
				if (x[i]._rN < 0 && x[i]._rN > -amax)
				{
					if (x[i]._rN > sol)
						sol = x[i]._rN;
				}
			}
		}
		
		Real ap2, t23, t45;
		ap2 = sol;
		t23 = (a0 * a0 - 2 * amax * amax + 2 * ap2 * ap2 - 2 * jmax*v0 + 2 * jmax*vt) / (2 * amax*jmax);
		t45 = (-a0 * a0 + 2 * a0*amax - 2 * amax * amax + 4 * amax*ap2 + 2 * jmax*tsyn*amax - 2 * ap2 * ap2 + 2 * jmax*v0 - 2 * jmax*vt) / (2 * amax*jmax);

		cout << "ap2 : " << ap2 << endl;
		cout << "t23 : " << t23 << endl;
		cout << "t45 : " << t45 << endl;
	}

};

//////////////////////// case 3 ////////////////////////
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

	void closedform()
	{
		Real t23, t67;
		t67 = -(2 * jmax*vt - 2 * jmax*v0 - 2 * a0*amax + a0 * a0 + 6 * amax *amax + 2 * sqrt((a0 * a0 * a0 * amax) / 3 -
			2 * a0*std::pow(amax,3) - (std::pow(a0, 4)) / 4 + std::pow(amax, 4) + a0 * a0 * amax * amax - jmax * jmax * v0 * v0 - jmax * jmax * vt * vt + 
			a0 * a0 * jmax*v0 - a0 * a0 * jmax*vt + 2 * jmax * jmax * v0*vt + amax * amax * jmax * jmax * tsyn * tsyn + 
			4 * amax*jmax*jmax * p0 - 4 * amax*jmax*jmax * pt - 2 * amax*amax*amax * jmax*tsyn - 2 * a0*amax*jmax*v0 +
			2 * a0*amax*jmax*vt + 2 * a0*amax *amax * jmax*tsyn - a0 *a0 * amax*jmax*tsyn + 2 * amax*jmax*jmax * tsyn*v0 +
			2 * amax*jmax*jmax * tsyn*vt)  - 2 * amax*jmax*tsyn) / (4 * amax*jmax);

		t23 = (a0 * a0 - 2 * jmax*v0 + 2 * jmax*vt + 2 * amax*jmax*t67) / (2 * amax*jmax);
		cout << t23 << endl;
		cout << t67 << endl;
	}
};

//////////////////////// case 4 ////////////////////////
class Step2PosTriZeroNegTrap // unknown : ap1, t34, t56
{
public:
	Real tsyn, a0, vt, v0, pt, p0, jmax, amax, vmax;

public:
	Step2PosTriZeroNegTrap(const Real _tsyn, const Real _a0, const Real _vt, const Real _v0, const Real _pt, const Real _p0, const Real _jmax, const Real _amax, const Real _vmax)
		: tsyn(_tsyn), a0(_a0), vt(_vt), v0(_v0), pt(_pt), p0(_p0), jmax(_jmax), amax(_amax), vmax(_vmax) {}

	void closedform()
	{
		Real a, b, c, d, e;
		a = 12;
		b = 24 * amax;
		c = -12 * a0 * a0 + 12 * amax * amax - 24 * jmax*vt - 24 * a0*amax - 24 * amax*jmax*tsyn + 24 * jmax*v0;
		d = 0;
		e = 4 * std::pow(a0, 3) * amax + 3 * std::pow(a0, 4) - 6 * a0 * a0 * amax * amax + 12 * jmax * jmax * v0 * v0 + 12 * jmax * jmax * vt * vt - 12 * a0 * a0 * jmax*v0 + 
			12 * amax * amax * jmax*v0 + 12 * a0 * a0 * jmax*vt - 12 * amax * amax * jmax*vt - 24 * jmax * jmax * v0*vt - 24 * amax*jmax * jmax * p0 + 
			24 * amax*jmax * jmax * pt + 12 * a0 * a0 * amax*jmax*tsyn - 24 * amax*jmax * jmax * tsyn*v0;

		ComplexNumber x[4];
		calcQuarticAlgebraicEqn(a, b, c, d, e, x);

		Real ep = 1E-8, sol = RealMax;
		for (int i = 0; i < 4; i++)
		{
			if (abs(x[i]._iN) < ep)
			{
				if (x[i]._rN > 0 && x[i]._rN < amax)
				{
					if (x[i]._rN < sol)
						sol = x[i]._rN;
				}
			}
		}

		//cout << "x[0] : " << x[0] << endl;
		//cout << "x[1] : " << x[1] << endl;
		//cout << "x[2] : " << x[2] << endl;
		//cout << "x[3] : " << x[3] << endl;

		Real ap1, t34, t56;
		ap1 = sol;
		t56 = -(a0 * a0 + 2 * amax * amax - 2 * ap1 * ap1 - 2 * jmax*v0 + 2 * jmax*vt) / (2 * amax*jmax);
		t34 = (a0 * a0 + 2 * a0*amax - 2 * amax * amax - 4 * amax*ap1 + 2 * jmax*tsyn*amax - 2 * ap1 * ap1 - 2 * jmax*v0 + 2 * jmax*vt) / (2 * amax*jmax);

		cout << "ap1 : " << ap1 << endl;
		cout << "t34 : " << t34 << endl;
		cout << "t56 : " << t56 << endl;
	}
};


//////////////////////// case 5 ////////////////////////
class Step2PosTrapZeroPosTrap // unknown : t23, t45, t67
{
public:
	Real tsyn, a0, vt, v0, pt, p0, jmax, amax, vmax;

public:
	Step2PosTrapZeroPosTrap(const Real _tsyn, const Real _a0, const Real _vt, const Real _v0, const Real _pt, const Real _p0, const Real _jmax, const Real _amax, const Real _vmax)
		: tsyn(_tsyn), a0(_a0), vt(_vt), v0(_v0), pt(_pt), p0(_p0), jmax(_jmax), amax(_amax), vmax(_vmax) {}

	void closedform()
	{
		Real a, b;
		a = -(24 * a0*std::pow(amax,3) - 8 * std::pow(a0, 3) * amax + 3 * std::pow(a0, 4) - 24 * std::pow(amax, 4) - 6 * a0 * a0 * amax * amax +
			12 * jmax * jmax * v0 * v0 + 12 * jmax * jmax * vt * vt - 12 * a0 * a0 * jmax*v0 + 12 * amax * amax * jmax*v0 +
			12 * a0 * a0 * jmax*vt - 12 * amax * amax * jmax*vt - 24 * jmax * jmax * v0*vt -
			24 * amax*jmax * jmax * p0 + 24 * amax*jmax * jmax * pt + 24 * std::pow(amax,3) * jmax*tsyn +
			24 * a0*amax*jmax*v0 - 24 * a0*amax*jmax*vt - 24 * amax*jmax * jmax * tsyn*vt);
		b = 24 * amax * amax * jmax * jmax * tsyn + 24 * amax*jmax * jmax * v0 -
			24 * std::pow(amax,3) * jmax + 24 * a0*amax * amax * jmax - 12 * a0 * a0 * amax*jmax - 24 * amax*jmax * jmax * vt;

		Real t23, t45, t67;
		t67 = a / b;
		t23 = -(-a0 * a0 + 4 * amax * amax + 2 * jmax*t67*amax + 2 * jmax*v0 - 2 * jmax*vt) / (2 * amax*jmax);
		t45 = (-a0 * a0 + 2 * a0*amax - 4 * amax * amax + 2 * jmax*tsyn*amax + 2 * jmax*v0 - 2 * jmax*vt) / (2 * amax*jmax);

		cout << "t23 : " << t23 << endl;
		cout << "t45 : " << t45 << endl;
		cout << "t67 : " << t67 << endl;
	}
};


//////////////////////// case 6 ////////////////////////
class Step2PosTrapZeroPosTri // unknown : t23, t45, ap1
{
public:
	Real tsyn, a0, vt, v0, pt, p0, jmax, amax, vmax;

public:
	Step2PosTrapZeroPosTri(const Real _tsyn, const Real _a0, const Real _vt, const Real _v0, const Real _pt, const Real _p0, const Real _jmax, const Real _amax, const Real _vmax)
		: tsyn(_tsyn), a0(_a0), vt(_vt), v0(_v0), pt(_pt), p0(_p0), jmax(_jmax), amax(_amax), vmax(_vmax) {}

	void closedform()
	{
		Real a, b, c, d, e;
		a = 12;
		b = -24 * amax;
		c = -12 * a0 * a0 - 12 * amax * amax + 24 * jmax*v0 - 24 * jmax*vt + 24 * amax*jmax*tsyn + 24 * a0*amax;
		d = 0;
		e = 3 * std::pow(a0, 4) - 8 * std::pow(a0, 3) * amax + 6 * a0 * a0 * amax * amax + 12 * jmax * jmax * v0 * v0 + 12 * jmax * jmax * vt * vt - 
			12 * a0 * a0 * jmax*v0 - 12 * amax * amax * jmax*v0 + 12 * a0 * a0 * jmax*vt + 12 * amax * amax * jmax*vt - 
			24 * jmax * jmax * v0*vt - 24 * amax*jmax * jmax * p0 + 24 * amax*jmax * jmax * pt + 24 * a0*amax*jmax*v0 - 
			24 * a0*amax*jmax*vt - 24 * amax*jmax * jmax * tsyn*vt;

		ComplexNumber x[4];
		calcQuarticAlgebraicEqn(a, b, c, d, e, x);

		Real ep = 1E-8, sol = RealMax;
		for (int i = 0; i < 4; i++)
		{
			if (abs(x[i]._iN) < ep)
			{
				if (x[i]._rN > 0 && x[i]._rN < amax)
				{
					if (x[i]._rN < sol)
						sol = x[i]._rN;
				}
			}
		}

		Real t23, t45, ap1;
		ap1 = sol;
		t23 = -(-a0 * a0 + 2 * amax * amax + 2 * ap1 * ap1 + 2 * jmax*v0 - 2 * jmax*vt) / (2 * amax*jmax);
		t45 = (-a0 * a0 + 2 * a0*amax - 2 * amax * amax - 4 * amax*ap1 + 2 * jmax*tsyn*amax + 2 * ap1 * ap1 + 2 * jmax*v0 - 2 * jmax*vt) / (2 * amax*jmax);
		
		cout << "ap1 : " << ap1 << endl;
		cout << "t23 : " << t23 << endl;
		cout << "t45 : " << t45 << endl;
	}
};


//////////////////////// case 7 ////////////////////////
class Step2PosTriZeroPosTrap // unknown : ap1, t34, t56
{
public:
	Real tsyn, a0, vt, v0, pt, p0, jmax, amax, vmax;

public:
	Step2PosTriZeroPosTrap(const Real _tsyn, const Real _a0, const Real _vt, const Real _v0, const Real _pt, const Real _p0, const Real _jmax, const Real _amax, const Real _vmax)
		: tsyn(_tsyn), a0(_a0), vt(_vt), v0(_v0), pt(_pt), p0(_p0), jmax(_jmax), amax(_amax), vmax(_vmax) {}

	void closedform()
	{
		Real a, b, c, d, e;
		a = 12;
		b = -24 * amax;
		c = -12 * a0 * a0 - 12 * amax * amax + 24 * jmax*v0 + 24 * a0*amax + 
			24 * amax*jmax*tsyn - 24 * jmax*vt;
		d = 0;
		e = 3 * std::pow(a0, 4) - 4 * std::pow(a0, 3) * amax + 6 * a0 * a0 * amax * amax + 12 * jmax * jmax * v0 * v0 + 
			12 * jmax * jmax * vt * vt - 12 * a0 * a0 * jmax*v0 - 12 * amax * amax * jmax*v0 + 
			12 * a0 * a0 * jmax*vt + 12 * amax * amax * jmax*vt - 24 * jmax * jmax * v0*vt + 24 * amax*jmax * jmax * p0 - 
			24 * amax*jmax * jmax * pt - 12 * a0 * a0 * amax*jmax*tsyn + 24 * amax*jmax * jmax * tsyn*v0;

		ComplexNumber x[4];
		calcQuarticAlgebraicEqn(a, b, c, d, e, x);

		Real ep = 1E-8, sol = RealMax;
		for (int i = 0; i < 4; i++)
		{
			if (abs(x[i]._iN) < ep)
			{
				if (x[i]._rN > 0 && x[i]._rN < amax)
				{
					if (x[i]._rN < sol)
						sol = x[i]._rN;
				}
			}
		}

		Real ap1, t34, t56;
		ap1 = sol;
		t34 = (-a0 * a0 + 2 * a0*amax - 2 * amax * amax - 4 * amax*ap1 + 2 * jmax*tsyn*amax + 2 * ap1 * ap1 + 2 * jmax*v0 - 2 * jmax*vt) / (2 * amax*jmax);
		t56 = -(-a0 * a0 + 2 * amax * amax + 2 * ap1 * ap1 + 2 * jmax*v0 - 2 * jmax*vt) / (2 * amax*jmax);

		cout << "ap1 : " << ap1 << endl;
		cout << "t34 : " << t34 << endl;
		cout << "t56 : " << t56 << endl;
	}
};


//////////////////////// case 8 ////////////////////////
class Step2PosTriZeroPosTri : public Function // unknown : ap1, ap2
{
public:
	Real tsyn, a0, vt, v0, pt, p0, jmax, amax, vmax;

public:
	Step2PosTriZeroPosTri(const Real _tsyn, const Real _a0, const Real _vt, const Real _v0, const Real _pt, const Real _p0, const Real _jmax, const Real _amax, const Real _vmax)
		: tsyn(_tsyn), a0(_a0), vt(_vt), v0(_v0), pt(_pt), p0(_p0), jmax(_jmax), amax(_amax), vmax(_vmax) {}

	VectorX func(const VectorX& x) const
	{
		VectorX val(2);
		Real ap1 = x(0), ap2 = x(1);
		val(0) = (-a0 * a0 + 2 * ap1 * ap1 + 2 * ap2 * ap2 + 2 * jmax*v0 - 2 * jmax*vt)/ (2 * jmax);
		val(1) = (6 * a0*ap1 * ap1 + 6 * jmax * jmax * p0 - 6 * jmax * jmax * pt - std::pow(a0, 3) - 6 * std::pow(ap1, 3) + 6 * std::pow(ap2, 3) + 
			6 * jmax * jmax * tsyn*v0 - 3 * a0 * a0 * jmax*tsyn + 6 * ap1 * ap1 * jmax*tsyn) / (6 * jmax * jmax);
		return val;
	}

	MatrixX InverseJacobian(const VectorX & x) const
	{
		MatrixX j(2, 2);
		Real ap1 = x(0), ap2 = x(1);
		Real det = -(2 * ap1*ap2*(2 * a0 - 3 * ap1 - 3 * ap2 + 2 * jmax*tsyn)) / std::pow(jmax, 3);
		j(0, 0) = (3 * ap2 * ap2) / (jmax * jmax) / det;
		j(1, 1) = (2 * ap1 / jmax) / det;
		j(0, 1) = -((2 * ap2) / jmax) / det;
		j(1, 0) = -((ap1*(2 * a0 - 3 * ap1 + 2 * jmax*tsyn)) / (jmax * jmax)) / det;
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
		vec2(0) = makeRandLU(0, 30); vec2(1) = makeRandLU(-30, 0);
		//cout << "initial ap1 : " << vec2(0) << ", initial ap2 : " << vec2(1) << endl;
		solver.setxN(2); solver.setfN(2); solver.setfunction(obj1);
		//time = clock();
		//for (int i = 0; i < 100; i++)
		solver.solve(vec2);
		//cout << "Computation time : " << clock() - time << endl;

		if (std::pow((10 - solver.getResultX()(0)), 2) > 1E-5)
		{
			errorcnt++;
			r = false;
		}
		cout << "initial ap1 : " << vec2(0) << ", initial ap2 : " << vec2(1) <<  "/ Result : " << solver.getResultX()(0) << ", " << solver.getResultX()(1) << ", " << r << endl;
	}

	Real success = (Real)(numOfExp - errorcnt) / (Real)(numOfExp)* 100;
	cout << "success per : " << success << "%" << endl;

	// PosTrapZeroNegTri Case Experiment --> closed-form
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
	obj2->closedform();
	cout << endl;

	// PosTrapZeroNegTrap Case Experiment --> closed-form
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
	obj3->closedform();
	cout << endl;

	// PosTriZeroNegTrap Case Experiment --> closed-form
	cout << "------- case 4 : PosTriZeroNegTrap -------" << endl;
	Step2PosTriZeroNegTrap obj4(6, 5, -28.75, 0, -27.7083, 0, 10, 15, 8.75);
	obj4.closedform();
	cout << endl;

	// PosTrapZeroPosTrap Case Experiment --> closed-form
	cout << "------- case 5 : PosTrapZeroPosTrap -------" << endl;
	Step2PosTrapZeroPosTrap obj5(7, 5, 58.75, 0, 216.0417, 0, 10, 15, 100);
	obj5.closedform();
	cout << endl;

	// PosTrapZeroPosTri Case Experiment --> closed-form
	cout << "------- case 6 : PosTrapZeroPosTri -------" << endl;
	Step2PosTrapZeroPosTri obj6(5.5, 5, 38.75, 0, 130.4167, 0, 10, 15, 100);
	obj6.closedform();
	cout << endl;

	// PosTriZeroPosTrap Case Experiment --> closed-form
	cout << "------- case 7 : PosTriZeroPosTrap -------" << endl;
	Step2PosTriZeroPosTrap obj7(5.5, 5, 38.75, 0, 95.4167, 0, 10, 15, 100);
	obj7.closedform();
	cout << endl;

	// PosTriZeroPosTri Case Experiment --> numerical
	cout << "------- case 8 : PosTriZeroPosTri -------" << endl;
	std::shared_ptr<Step2PosTriZeroPosTri> obj8 = std::shared_ptr<Step2PosTriZeroPosTri>(new Step2PosTriZeroPosTri(4, 5, 18.75, 0, 39.7917, 0, 10, 15, 100));
	vec2(0) = makeRandLU(0, 15); vec2(1) = makeRandLU(0, 15);
	cout << "initial ap1 : " << vec2(0) << endl;
	cout << "initial ap2 : " << vec2(1) << endl;
	solver.setxN(2); solver.setfN(2); solver.setfunction(obj8);
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