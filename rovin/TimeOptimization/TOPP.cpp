#include "TOPP.h"

using namespace std;

namespace rovin {

	TOPP::TOPP(const MatrixX & q_data, const SerialOpenChainPtr & soc, 
		const Real ds, const Real vi, const Real vf, const Real si, const Real sf)
	{
		// _q_data dimesion : dof * number of data
		_q_data = q_data;
		_soc = soc;
		_ds = ds;
		_dof = _soc->getNumOfJoint();
		_vi = vi;
		_vf = vf;
		_si = si;
		_sf = sf;

		// insert torque constraint, velocity constraint
		_torqueConstraint.resize(_dof * 2);
		_velConstraint.resize(_dof * 2);
		_accConstraint.resize(_dof * 2);
		for (int i = 0; i < _dof; i++)
		{
			_torqueConstraint[i] = _soc->getMotorJointPtr(i)->getLimitTorqueUpper();
			_torqueConstraint[i + _dof] = _soc->getMotorJointPtr(i)->getLimitTorqueLower();

			_velConstraint[i] = _soc->getMotorJointPtr(i)->getLimitVelUpper();
			_velConstraint[i + _dof] = _soc->getMotorJointPtr(i)->getLimitVelLower();

			_accConstraint[i] = _soc->getMotorJointPtr(i)->getLimitAccUpper();
			_accConstraint[i] = _soc->getMotorJointPtr(i)->getLimitAccLower();
		}

		//_torqueConstraint *= 2.0;
		//cout << _torqueConstraint << endl;

		// make b-spline by wooyoung
		unsigned int degreeOfBSpline = 3;
		unsigned int orderOfBSpline = degreeOfBSpline + 1;
		unsigned int MaxNumOfCP = 7;
		unsigned int numData = _q_data.cols();
		if (numData > MaxNumOfCP)
		{
			cout << "Make b-spline using fitting" << endl;
			_q = BSplineFitting(_q_data, orderOfBSpline, MaxNumOfCP, si, sf);
		}
		else
		{
			cout << "Make b-spline using interpolation" << endl;
			_q = BSplineInterpolation(_q_data, orderOfBSpline, si, sf);
		}
		_dqds = _q.derivative();
		_ddqdds = _dqds.derivative();

		// B-spline experiment 1
		//VectorX knot;
		//unsigned int degreeOfBSpline = 4;
		//unsigned int numOfCP = 7;
		//MatrixX data(_dof, 7);
		//knot.resize(degreeOfBSpline + numOfCP + 1);
		//for (unsigned int i = 0; i < degreeOfBSpline + 1; i++)
		//{
		//	knot[i] = _si;
		//	knot[knot.size() - i - 1] = _sf;
		//}
		//Real delta = _sf / (numOfCP - degreeOfBSpline);
		//for (unsigned int i = 0; i < numOfCP - degreeOfBSpline; i++)
		//{
		//	knot[degreeOfBSpline + 1 + i] = delta * (i + 1);
		//}

		//_q = BSpline<-1, -1, -1>(knot, _q_data);
		//_dqds = _q.derivative();
		//_ddqdds = _dqds.derivative();

	}

	bool TOPP::checkMVCCondition(Real alpha, Real beta)
	{
		if (beta > alpha)
			return true;
		else
			return false;
	}

	VectorX TOPP::calculateA(Real s)
	{
		VectorX a(_dof * 2);
		VectorX q = _q(s);
		//cout << "q : " << q << endl;
		VectorX qs = _dqds(s);
		//cout << "qs : " << qs << endl;
		VectorX qdot(q.size());
		qdot.setZero();

		StatePtr state = _soc->makeState();
		state->setJointStatePos(q);
		state->setJointStateVel(qdot);
		state->setJointStateAcc(qs);

		_soc->solveInverseDynamics(*state);

		VectorX tau = state->getJointStateTorque();
		//cout << "torque : " << tau << endl;

		qs.setZero();
		state->setJointStatePos(q);
		state->setJointStateVel(qdot);
		state->setJointStateAcc(qs);
		_soc->solveInverseDynamics(*state);

		VectorX tmp_a = tau - state->getJointStateTorque();
		//cout << "tmp torque : " << tmp_a << endl;
		
		for (int i = 0; i < _dof; i++)
		{
			a[i] = tmp_a[i];
			a[i + _dof] = -tmp_a[i];
		}

		//cout << "a : " << a << endl;
		return a;
	}

	std::vector<VectorX> TOPP::calculateBandC(Real s)
	{
		std::vector<VectorX> result;
		VectorX b(_dof * 2);
		VectorX c(_dof * 2);

		VectorX q = _q(s);
		VectorX qdot = _dqds(s);
		VectorX qddot = _ddqdds(s);
		StatePtr state = _soc->makeState();

		state->setJointStatePos(q);
		state->setJointStateVel(qdot);
		state->setJointStateAcc(qddot);

		_soc->solveInverseDynamics(*state);

		VectorX tau = state->getJointStateTorque();

		state->setJointStatePos(q);
		state->setJointStateVel(qdot.setZero());
		state->setJointStateAcc(qddot.setZero());

		_soc->solveInverseDynamics(*state);

		VectorX tmp_c = state->getJointStateTorque();
		VectorX tmp_b = tau - tmp_c;
		
		for (int i = 0; i < _dof; i++)
		{
			//b[i] = -tmp_b[i];
			//b[i + _dof] = tmp_b[i];
			//c[i] = -tmp_c[i];
			//c[i + _dof] = tmp_c[i];
			b[i] = tmp_b[i];
			b[i + _dof] = -tmp_b[i];
			c[i] = tmp_c[i];
			c[i + _dof] = -tmp_c[i];
		}

		for (int i = 0; i < _dof; i++)
		{
			//c[i] += _torqueConstraint[i];
			//c[i + _dof] -= _torqueConstraint[i + _dof];
			c[i] -= _torqueConstraint[i];
			c[i + _dof] += _torqueConstraint[i + _dof];
		}
		
		result.push_back(b);
		result.push_back(c);

		//cout << "b : " << b << endl;
		//cout << "c : " << c << endl;
		//cout << "b+c : " << b + c << endl;

		return result;
	}

	Vector2 TOPP::determineAlphaBeta(Real s, Real sdot)
	{
		VectorX a = calculateA(s);
		//cout << "a : " << a << endl;

		// zero-inertia point check
		bool zero_inertia_swi = false;
		for (int i = 0; i < a.size(); i++)
		{
			//cout << a[i] << endl;
			if (std::abs(a[i]) < RealEps)
			{
				zero_inertia_swi = true;
				break;
			}
		}

		if (true)//!zero_inertia_swi)
		{
			VectorX q = _q(s);
			VectorX qdot = _dqds(s)*sdot;
			VectorX qddot = _ddqdds(s)*sdot*sdot;
			StatePtr state = _soc->makeState();

			state->setJointStatePos(q);
			state->setJointStateVel(qdot);
			state->setJointStateAcc(qddot);

			_soc->solveInverseDynamics(*state);

			VectorX tau = state->getJointStateTorque();
			//cout << "tau : " << tau << endl;
			VectorX left_vec(_dof * 2);
			for (int i = 0; i < _dof; i++)
			{
				left_vec(i) = -tau(i) + _torqueConstraint(i);
				left_vec(i + _dof) = tau(i) - _torqueConstraint(i + _dof);
				//a_agg[i] = a[i];
				//a_agg[i + _dof] = a[i];
			}

			//cout << "a : " << a << endl;
			//cout << "left_vec(b*sdot^(2) + c) : " << left_vec << endl;

			Vector2 result; // result[0] is alpha, result[1] is beta
			result[0] = -std::numeric_limits<Real>::max();
			result[1] = std::numeric_limits<Real>::max();

			for (int i = 0; i < _dof * 2; i++)
			{
				Real tmp = left_vec(i) / a(i);
				if (a[i] > RealEps) // upper bound beta
				{
					if (tmp < result[1])
						result[1] = tmp;
				}
				else if (a[i] <-RealEps)// lower bound alpha
				{
					if (tmp > result[0])
						result[0] = tmp;
				}
		//else // lower bound alpha
		//{
		//	if (tmp > result[0])
		//		result[0] = tmp;
		//}
			}

			//cout << "result : " << result << endl;

			return result;
		}
		else // Dynamic  singularity
		{
			// TODO
			return Vector2();
		}
	}

	Real TOPP::calculateMVCPoint(Real s)
	{
		Real sdot_MVC = std::numeric_limits<Real>::max();

		Vector2 alphabeta = determineAlphaBeta(s, 0);
		if (alphabeta[0] > alphabeta[1])
			return 0;

		VectorX a = calculateA(s);
		//cout << endl;
		//cout << "a : " << a << endl;

		std::vector<VectorX> BandC = calculateBandC(s);
		VectorX b = BandC[0];
		VectorX c = BandC[1];

		
		LOGIF(((a.size() == b.size()) && (a.size() == c.size()) && (b.size() == c.size())), "TOPP::calulateMVCPoint error : a, b, c vector size is wrong.");
		int nconstraints = _dof * 2;

		unsigned int kk;
		unsigned int mm;

		for (int k = 0; k < nconstraints; k++)
		{
			for (int m = k + 1; m < nconstraints; m++)
			{
				Real num, denum, r;

				// If we have a pair of alpha and beta bounds, then determine the sdot for which the two bounds are equal
				if (a(k)*a(m) < 0)
				{
					num = a(k) * c(m) - a(m) * c(k);
					denum = -a(k) * b(m) + a(m) * b(k);
					if (std::abs(denum) > 1e-10)
					{
						r = num / denum;
						if (r >= 0)
						{
							sdot_MVC = std::min(sdot_MVC, sqrt(r));
							kk = k;
							mm = m;
						}
							
					}
				}
			}
		}
		return sdot_MVC;
	}

	void TOPP::calculateFinalTime()
	{
		int integrateType = 0; // 0 : Wooyoung, 1 : GQ, 2 : Euler, 3 : Trapz

		if (integrateType == 0)
		{
			std::list<Real>::iterator it_sdot = ++(_sdot.begin());

			Real s_k = _s.front(), sdot_k = _sdot.front();
			Real s_k1, sdot_k1;
			Real sum = 0;
			for (std::list<Real>::iterator it_s = ++(_s.begin()); it_s != (_s.end()); ++it_s)
			{
				s_k1 = *(it_s); sdot_k1 = *(it_sdot);
				sum += 2 * (s_k1 - s_k) / (sdot_k1 + sdot_k);
				s_k = s_k1; sdot_k = sdot_k1;
				it_sdot++;
			}
			_tf_result = sum;
		}
		else if (integrateType == 1)
		{
			int numOfGQPoint = 40;
			GaussianQuadrature GQ(numOfGQPoint, _si, _sf);

			MatrixX reverse_sdot;
			if (_sdot.front() < RealEps)
			{
				reverse_sdot = MatrixX(1, _sdot.size() - 2);
				_sdot.pop_front();
				_s.pop_front();
			}
			else
				reverse_sdot = MatrixX(1, _sdot.size() - 1);

			int cnt = 0;
			for (std::list<Real>::iterator it = (_sdot.begin()); it != --(_sdot.end()); ++it)
			{
				reverse_sdot(0, cnt) = 1.0 / (*it);
				cnt++;
			}
		
			unsigned int degreeOfBSpline = 3;
			unsigned int orderOfBSpline = degreeOfBSpline + 1;
			unsigned int MaxNumOfCP = 10;
			BSpline<-1, -1, -1> f;
			f = BSplineFitting(reverse_sdot, orderOfBSpline, MaxNumOfCP, _s.front(), _s.back());

			VectorX sum = VectorX::Zero(1);
			const VectorX& weights = GQ.getWeights();
			const VectorX& QueryPoints = GQ.getQueryPoints();

			for (int i = 0; i < numOfGQPoint; i++)
			{
				sum(0) += f(QueryPoints(i))(0,0) * weights(i);
			}
			_tf_result = sum(0);
		}
		else if (integrateType == 2)
		{
			Real sum = 0;
			Real s_cur;
			if (_sdot.front() < RealEps)
			{
				_sdot.pop_front();
				_s.pop_front();
			}
			while (!_s.empty())
			{
				s_cur = _s.front();
				_s.pop_front();
				if (!_s.empty())
					sum += 1 / (_sdot.front()) * (_s.front() - s_cur);
				_sdot.front();
			}
			_tf_result = sum;
		}
		else if (integrateType == 3)
		{
			if (_sdot.front() < RealEps)
			{
				_sdot.pop_front();
				_s.pop_front();
			}
			Real sum = 0;;
			std::vector<Real> reverse_sdot;
			std::vector<Real> s;
			std::list<Real>::iterator it_s = _s.begin();
			for (std::list<Real>::iterator it = (_sdot.begin()); it != --(_sdot.end()); ++it)
			{
				reverse_sdot.push_back(1.0 / (*it));
				s.push_back((*it_s));
				it_s++;
			}

			for (int i = 0; i < reverse_sdot.size() - 1; i++)
			{
				sum += 0.5 * (reverse_sdot[i + 1] + reverse_sdot[i]) * (s[i + 1] - s[i]);
			}
			_tf_result = sum;
		}
	}

	void TOPP::calculateTorqueTrajectory()
	{
		// calculate time and joint angle
		VectorX t(_s.size());
		t(0) = 0;

		std::list<Real>::iterator it_s = ++(_s.begin());
		std::list<Real>::iterator it_sdot = ++(_sdot.begin());
		Real s_k = _s.front(), sdot_k = _sdot.front();
		Real s_k1, sdot_k1, dt;

		for (int i = 1; i < _s.size(); i++)
		{
			s_k1 = *(it_s); sdot_k1 = *(it_sdot);
			dt = 2 * (s_k1 - s_k) / (sdot_k1 + sdot_k);
			t(i) = t(i - 1) + dt;
			s_k = s_k1; sdot_k = sdot_k1;
			it_s++; it_sdot++;
		}

		MatrixX q_data(_dof, t.size());
		it_s = _s.begin();
		for (int i = 0; i < t.size(); i++)
			q_data.col(i) = _q((*it_s)++);

		// make b-spline
		// TODO --> Woo young
		BSpline<-1, -1, -1> q;
		BSpline<-1, -1, -1> qdot = q.derivative();
		BSpline<-1, -1, -1> qddot = qdot.derivative();

		// solve inverse dynamics
		_torque_result = MatrixX(_dof, t.size());
		StatePtr state = _soc->makeState();
		for (int i = 0; i < t.size(); i++)
		{
			state->setJointStatePos(q(t(i)));
			state->setJointStateVel(qdot(t(i)));
			state->setJointStateAcc(qddot(t(i)));
			_soc->solveInverseDynamics(*state);
			_torque_result.col(i) = state->getJointStateTorque();
		}
	}

	Real TOPP::calculateMVCPointExclude(Real s, int iExclude)
	{
		Real sdot_MVC = std::numeric_limits<Real>::max();

		Vector2 alphabeta = determineAlphaBeta(s, 0);
		if (alphabeta[0] > alphabeta[1])
			return 0;

		VectorX a = calculateA(s);
		std::vector<VectorX> BandC = calculateBandC(s);
		VectorX b = BandC[0];
		VectorX c = BandC[1];

		LOGIF(((a.size() == b.size()) && (a.size() == c.size()) && (b.size() == c.size())), "TOPP::calulateMVCPoint error : a, b, c vector size is wrong.")
		int nconstraints = _dof * 2;

		for (int k = 0; k < nconstraints; k++)
		{
			for (int m = k + 1; m < nconstraints; m++)
			{
				// exclude iExclude-th inequality constraint
				if (k == iExclude || m == iExclude) {
					continue;
				}

				Real num, denum, r;

				// If we have a pair of alpha and beta bounds, then determine the sdot for which the two bounds are equal
				if (a(k)*a(m) < 0)
				{
					num = a(k) * c(m) - a(m) * c(k);
					denum = -a(k) * b(m) + a(m) * b(k);
					if (std::abs(denum) > 1e-10)
					{
						r = num / denum;
						if (r >= 0)
							sdot_MVC = std::min(sdot_MVC, sqrt(r));
					}
				}
			}
		}

		return sdot_MVC;
	}

	void TOPP::farwardIntegrate(Real & s, Real & sdot, Real sddot)
	{		
		// Explicit Euler
		//Real tmp = 2 * sddot*_ds + sdot*sdot;
		//LOGIF((tmp > 0), "TOPP::farwardIntegrate error : the value has to be positive.");
		//s = s + _ds;
		//sdot = sqrt(tmp);

		// prediction-correction method
		/* step 1 : prediction using explicit euler */
		//Real s_p = s + _ds;
		//Real tmp = 2 * sddot*_ds + sdot*sdot;
		//LOGIF((tmp > 0), "TOPP::farwardIntegrate error : the value has to be positive.");
		//Real sdot_p = sqrt(tmp);
		//Real sddot_p = determineAlphaBeta(s_p, sdot_p)(1);
		///* step 2 : prediction using trapz */
		//tmp = (sddot + sddot_p)*_ds + sdot*sdot;
		//LOGIF((tmp > 0), "TOPP::farwardIntegrate error : the value has to be positive.");
		//s = s + _ds;
		//sdot = sqrt(tmp);

		// RK4
		Real s_RK = s + _ds;
		Real s_RK_half = s + 0.5*_ds;
		Real f = 2 * sddot;
		Real tmp = sdot*sdot + 0.5 * _ds * f;
		LOGIF((tmp > 0), "TOPP::farwardIntegrate error : the value has to be positive.");
		Real sdot_RK_1 = sqrt(tmp);
		Real f1 = 2 * determineAlphaBeta(s_RK_half, sdot_RK_1)(1);
		tmp = sdot*sdot + 0.5 * _ds * f1;
		LOGIF((tmp > 0), "TOPP::farwardIntegrate error : the value has to be positive.");
		Real sdot_RK_2 = sqrt(tmp);
		Real f2 = 2 * determineAlphaBeta(s_RK_half, sdot_RK_2)(1);
		tmp = sdot*sdot + _ds * f2;
		LOGIF((tmp > 0), "TOPP::farwardIntegrate error : the value has to be positive.");
		Real sdot_RK_3 = sqrt(tmp);
		tmp = sdot*sdot + _ds * (1.0/6.0*sddot + 1.0/3.0*f1 + 1.0/3.0*f2 + 1.0/6.0*2* determineAlphaBeta(s_RK, sdot_RK_3)(1));
		LOGIF((tmp > 0), "TOPP::farwardIntegrate error : the value has to be positive.");
		sdot = sqrt(tmp);
		s = s + _ds;
	}

	void TOPP::backwardIntegrate(Real & s, Real & sdot, Real sddot)
	{
		// Explicit Euler
		//Real tmp = -2 * sddot*_ds + sdot*sdot;
		//LOGIF((tmp > 0), "TOPP::backwardIntegrate error : the value has to be positive.");
		//s = s - _ds;
		//sdot = sqrt(tmp);

		// prediction-correction method
		/* step 1 : prediction using explicit euler */
		//Real s_p = s - _ds;
		//Real tmp = - 2 * sddot*_ds + sdot*sdot;
		//LOGIF((tmp > 0), "TOPP::backwardIntegrate error : the value has to be positive.");
		//Real sdot_p = sqrt(tmp);
		//Real sddot_p = determineAlphaBeta(s_p, sdot_p)(0);
		///* step 2 : prediction using trapz */
		//tmp = -(sddot + sddot_p)*_ds + sdot*sdot;
		//LOGIF((tmp > 0), "TOPP::backwardIntegrate error : the value has to be positive.");
		//s = s - _ds;
		//sdot = sqrt(tmp);

		// RK4
		Real s_RK = s - _ds;
		Real s_RK_half = s - 0.5*_ds;
		Real f = 2 * sddot;
		Real tmp = sdot*sdot - 0.5 * _ds * f;
		LOGIF((tmp > 0), "TOPP::farwardIntegrate error : the value has to be positive.");
		Real sdot_RK_1 = sqrt(tmp);
		Real f1 = 2 * determineAlphaBeta(s_RK_half, sdot_RK_1)(0);
		tmp = sdot*sdot - 0.5 * _ds * f1;
		LOGIF((tmp > 0), "TOPP::farwardIntegrate error : the value has to be positive.");
		Real sdot_RK_2 = sqrt(tmp);
		Real f2 = 2 * determineAlphaBeta(s_RK_half, sdot_RK_2)(0);
		tmp = sdot*sdot - _ds * f2;
		LOGIF((tmp > 0), "TOPP::farwardIntegrate error : the value has to be positive.");
		Real sdot_RK_3 = sqrt(tmp);
		tmp = sdot*sdot - _ds * (1.0 / 6.0*sddot + 1.0 / 3.0*f1 + 1.0 / 3.0*f2 + 1.0 / 6.0 * 2 * determineAlphaBeta(s_RK, sdot_RK_3)(1));
		LOGIF((tmp > 0), "TOPP::farwardIntegrate error : the value has to be positive.");
		sdot = sqrt(tmp);
		s = s - _ds;
	}

	bool TOPP::findNearestSwitchPoint(Real s)
	{

		Real ds = 0.0005;
		
		Real s_bef = s;
		Real sdot_bef = calculateMVCPoint(s_bef);
		Real s_cur = s_bef + ds;
		Real sdot_cur = calculateMVCPoint(s_cur);
		Real s_next = s_bef + 2 * ds;
		Real sdot_next = calculateMVCPoint(s_next);

		// end criterion
		if (s_next >= _sf)
			return false;

		// ���� topp �ڵ忡�� tan_bef/tan_cur �� discontinous ����Ʈ �߰��ϴµ� �װ� ���� �Ȱ�
		// �ȳ־ �ɰͰ����� �׷��Ÿ� tan_bef/tan_cur ���� ���� �ʿ� ����
		//Real tan_bef = (sdot_cur - sdot_bef) / ds;
		//Real tan_cur = (sdot_next - sdot_cur) / ds;
		Real diff_bef = determineAlphaBeta(s_bef, sdot_bef)[0] / sdot_bef - (sdot_cur - sdot_bef) / ds;
		Real diff_cur = determineAlphaBeta(s_cur, sdot_cur)[0] / sdot_cur - (sdot_next - sdot_cur) / ds;

		VectorX a_bef = calculateA(s_bef);
		std::vector<VectorX> bc_bef = calculateBandC(s_bef);
		VectorX a_cur = calculateA(s_cur);
		std::vector<VectorX> bc_cur = calculateBandC(s_cur);


		bool singularFound = false;
		Real s_sing, a_sing, b_sing, c_sing, sdot_star, sdot_plus, s_save; ///> variables for finding singular points
		Real sdot_min = std::numeric_limits<Real>::max();
		int idx_save;

		while (true)
		{
			
			// step 1: singular point check
			// include calculating \lambda
			for (unsigned int i = 0; i < 2 * _dof; i++) // number of inequality
			{
				if (a_bef[i] * a_cur[i] <= 0)
				{
					// �� 0�� �ƴϸ� ���� �ڵ�ó�� �������� �ؼ� �� ���� s�� ã�� ������ �ʿ��ҵ�..
					// ���� ���� ���Ʒ� step 2, 3���� ����� �˰��� �ʿ��ҵ�
					Real adiff = a_cur[i] - a_bef[i];
					if (std::abs(adiff) > 1E-10)
					{
						// calculate '_sing' variables as internal dividing points
						// abc�� ������������?? ���� ���� s_sing ���� ���ϴ°� �ƴϰ�??
						Real tmpalpha = -a_bef[i] / adiff;
						s_sing = tmpalpha*s_cur + (1 - tmpalpha)*s_bef;
						a_sing = calculateA(s_sing)[i];
						std::vector<VectorX> bcTmp = calculateBandC(s_sing);
						b_sing = bcTmp[0](i);
						c_sing = bcTmp[1](i);
					}
					else
					{
						s_sing = s_cur;
						a_sing = a_cur[i];
						b_sing = bc_cur[0][i];
						c_sing = bc_cur[1][i];
					}

					if (b_sing > 0 && c_sing < 0)
					{
						sdot_star = sqrt(-c_sing / b_sing);
						sdot_plus = calculateMVCPointExclude(s_sing, i);

						if (sdot_star < sdot_plus && sdot_star < sdot_min)
						{
							singularFound = true;
							s_save = s_sing;
							sdot_min = sdot_star;
							idx_save = i;
						}

					}
				}
			}
			if (singularFound)
			{
				Real diffeps = 1E-5;
				Real ap, bp, cp; ///> p for prime, differentiated to s
				ap = (calculateA(s_save + diffeps)(idx_save) - calculateA(s_save)(idx_save)) / diffeps;
				std::vector<VectorX> bcTmp = calculateBandC(s_save);
				std::vector<VectorX> bcTmpEps = calculateBandC(s_save + diffeps);
				bp = (bcTmpEps[0](idx_save) - bcTmp[0](idx_save)) / diffeps;
				cp = (bcTmpEps[1](idx_save) - bcTmp[1](idx_save)) / diffeps;

				Real lambda;

				if ((2 * bcTmp[0](idx_save) + ap)*sdot_min > 1E-10)
					lambda = (-bp * sdot_min*sdot_min - cp) / ((2 * bcTmp[0](idx_save) + ap)*sdot_min);
				else
					lambda = 0.0; // copy of original code....

				SwitchPoint sw(s_save, sdot_min, SwitchPoint::SINGULAR, lambda);
				_switchPoint.push_back(sw);
				return true;
			}
			// step 1 -END-

			// step 2: tanget point check
			if ((diff_bef*diff_cur < 0) && (std::abs(diff_cur) < 1))
			{
				SwitchPoint sw(s_cur, sdot_cur, SwitchPoint::TANGENT, 0.0);
				_switchPoint.push_back(sw);
				return true;
			}
			// step 2 -END-


			// step 3: discontinuous point check
			if (std::abs(sdot_next - sdot_cur) > 1000 * std::abs(sdot_cur - sdot_bef))
			{
				if (sdot_next > sdot_cur)
				{
					SwitchPoint sw(s_cur, sdot_cur, SwitchPoint::DISCONTIUOUS, 0.0);
					_switchPoint.push_back(sw);
					return true;
				}
				else
				{
					SwitchPoint sw(s_next, sdot_next, SwitchPoint::DISCONTIUOUS, 0.0);
					_switchPoint.push_back(sw);
					return true;
				}
			}
			// step 3 -END-

			// if any kinds of switch point is not detected, proceed to next s
			s_bef = s_cur;
			sdot_bef = sdot_cur;
			s_cur = s_next;
			sdot_cur = sdot_next;
			s_next += ds;
			sdot_next = calculateMVCPoint(s_next);

			//tan_bef = tan_cur;
			//tan_cur = (sdot_next - sdot_cur) / ds;
			diff_bef = diff_cur;
			diff_cur = determineAlphaBeta(s_cur, sdot_cur)[0] / sdot_cur - (sdot_next - sdot_cur) / ds;

			a_bef = a_cur;
			bc_bef = bc_cur;
			a_cur = calculateA(s_cur);
			bc_cur = calculateBandC(s_cur);

			// end criterion
			if (s_next >= _sf)
				return false;
		}
	}

	void TOPP::generateTrajectory()
	{
		// initialization
		Real s_cur = 0;
		Real sdot_cur = _vi / _dqds(0.0001).norm();
		_s.push_back(s_cur);
		_sdot.push_back(sdot_cur);

		Vector2 alphabeta = determineAlphaBeta(s_cur, sdot_cur);
		Real alpha_cur = alphabeta(0);
		Real beta_cur = alphabeta(1);
		//cout << "alpha initial : " << alpha_cur << endl;
		//cout << "beta initial : " << beta_cur << endl;

		// list used when backward integration
		std::list<Real> _s_tmp;
		std::list<Real> _sdot_tmp;

		bool FI_SW = true; ///< forward integration switch
		bool BI_SW = false; ///< backward integration switch
		bool I_SW = true; ///< integration switch

		bool swiPoint_swi = false;
		unsigned int numOfSPInt = 1;
		
		while (I_SW)
		{
			//cout << "I_SW" << endl;
			// Forward integration
			while (FI_SW)
			{
				s_FI_jk.push_back(s_cur);
				sd_FI_jk.push_back(sdot_cur);

				//std::cout << "FI_SW" << endl;
				//std::cout << "s_cur : " << s_cur << endl;
				//std::cout << "sdot_cur : " << sdot_cur << endl;

				// forward intergration
				farwardIntegrate(s_cur, sdot_cur, beta_cur); ///< update s_cur, sdot_cur

				// save trajectory points
				_s.push_back(s_cur);
				_sdot.push_back(sdot_cur);

				// calculate alpha and beta
				alphabeta = determineAlphaBeta(s_cur, sdot_cur);
				alpha_cur = alphabeta(0);
				beta_cur = alphabeta(1);
				//std::cout << "alpha_cur : " << alphabeta(0) << endl;
				//std::cout << "beta_cur : " << alphabeta(1) << endl;

				//if (s_cur > 0.183)
				//{
				//	int aaaa = 3;
				//}

				//cout << alpha_cur << '\t' << beta_cur << endl;

				if (!checkMVCCondition(alpha_cur, beta_cur)) // case (a)
				{
					//saveRealVector2txt(s_FI_jk, "C:/Users/crazy/Desktop/Time optimization/s_sw.txt");
					//saveRealVector2txt(sd_FI_jk, "C:/Users/crazy/Desktop/Time optimization/sdot_sw.txt");

					// fine nearest switch point
					cout << "s_cur : " << s_cur << endl;
					swiPoint_swi = findNearestSwitchPoint(s_cur);

					cout << "switching point" << endl;
					cout << "switching point switch : " << swiPoint_swi << endl;
					cout << "switching point size : " << _switchPoint.size() << endl;
					if (_switchPoint.size() > 0)
					{
						for (int i = 0; i < _switchPoint.size(); i++)
						{
							cout << "switch point s : " << _switchPoint[i]._s << endl;
							cout << "switch point sdot : " << _switchPoint[i]._sdot << endl;
						}
						//cout << "switching point value : " << _switchPoint[_switchPoint.size() - 1]._s << ", " << _switchPoint[_switchPoint.size() - 1]._sdot << endl;
					}
					// if swtich point doesn't exist until s_end --> go to step 3
					if (!swiPoint_swi)
					{
						//std::cout << "Can not find switching point." << std::endl;
						FI_SW = false;
						BI_SW = false;
						I_SW = false;
					}
					else // if switch point exist --> go to backward integration
					{

						s_cur = _switchPoint[_switchPoint.size() - 1]._s;
						sdot_cur = _switchPoint[_switchPoint.size() - 1]._sdot;

						_s_tmp.push_front(s_cur);
						_sdot_tmp.push_front(sdot_cur);

						if (_switchPoint[_switchPoint.size() - 1]._id == SwitchPoint::SINGULAR)
						{
							Real lambda = _switchPoint[_switchPoint.size() - 1]._lambda;
							for (int i = 0; i < numOfSPInt; i++)
							{
								s_cur -= _ds;
								sdot_cur -= -lambda * _ds / s_cur;
								_s_tmp.push_front(s_cur);
								_sdot_tmp.push_front(sdot_cur);
							}
						}
						FI_SW = false;
						BI_SW = true;
					}
				}
				else if (s_cur > _sf) ///< go to step 3, case (a), s_cur �� ������ s_end �Ѿ�� ��!!
				{
					FI_SW = false;
					BI_SW = false;
					I_SW = false;
				}
				else if (sdot_cur < 1e-10) ///< Debugging ���, case (a)
				{
					FI_SW = false;
					BI_SW = false;
					I_SW = false;
				}
			}

			//saveRealVector2txt(s_FI_jk, "D:/jkkim/Documents/matlabTest/sFI.txt");
			//saveRealVector2txt(sd_FI_jk, "D:/jkkim/Documents/matlabTest/sdotFI.txt");


			// Backward intergration
			while (BI_SW)
			{
				s_BI_jk.push_back(s_cur);
				sd_BI_jk.push_back(sdot_cur);

				//std::cout << "BI_SW" << endl;
				//std::cout << "s_cur : " << s_cur << endl;
				//std::cout << "sdot_cur : " << sdot_cur << endl;

				// calculate alpha and beta
				alphabeta = determineAlphaBeta(s_cur, sdot_cur);
				alpha_cur = alphabeta(0);
				beta_cur = alphabeta(1);

				//cout << "alpha_cur : " << alpha_cur << endl;
				//cout << "beta_cur : " << beta_cur << endl;

				// backward integration
				backwardIntegrate(s_cur, sdot_cur, alpha_cur); ///< update s_cur, sdot_cur
				

				// save trajectory points
				_s_tmp.push_front(s_cur);
				_sdot_tmp.push_front(sdot_cur);

				if (s_cur <= _s.back())
				{
					//saveRealVector2txt(s_BI_jk, "C:/Users/crazy/Desktop/Time optimization/s_bsw.txt");
					//saveRealVector2txt(sd_BI_jk, "C:/Users/crazy/Desktop/Time optimization/sdot_bsw.txt");

					cout << "Enter if routine" << endl;

					_s_tmp.pop_front();
					_sdot_tmp.pop_front();

					sdot_cur = (sdot_cur - _sdot_tmp.front()) / (s_cur - _s_tmp.front()) * (_s.back() - s_cur) + sdot_cur;
					s_cur = _s.back();

					LOGIF(_sdot.back() > sdot_cur,"Backward intergration error:_sdot.back() has to be larger than sdot_cur.");

					_s_tmp.push_front(s_cur);
					_sdot_tmp.push_front(sdot_cur);

					while (_sdot.back() > sdot_cur)
					{
						s_BI_jk.push_back(s_cur);
						sd_BI_jk.push_back(sdot_cur);

						_s.pop_back();
						_sdot.pop_back();
						alphabeta = determineAlphaBeta(s_cur, sdot_cur);
						alpha_cur = alphabeta(0);

						//cout << "s_cur : " << s_cur << endl;
						//cout << "sdot_cur : " << sdot_cur << endl;
						//cout << "alpha_cur : " << alpha_cur << endl;

						backwardIntegrate(s_cur, sdot_cur, alpha_cur); ///< update s_cur, sdot_cur
						_s_tmp.push_front(s_cur);
						_sdot_tmp.push_front(sdot_cur);
					}

					s_BI_jk.push_back(s_cur);
					sd_BI_jk.push_back(sdot_cur);
					//saveRealVector2txt(s_BI_jk, "C:/Users/crazy/Desktop/Time optimization/s_bsw.txt");
					//saveRealVector2txt(sd_BI_jk, "C:/Users/crazy/Desktop/Time optimization/sdot_bsw.txt");

					_s_tmp.pop_front();
					_sdot_tmp.pop_front();

					_s.insert(_s.end(), _s_tmp.begin(), _s_tmp.end());
					_sdot.insert(_sdot.end(), _sdot_tmp.begin(), _sdot_tmp.end());

					_s_tmp.clear();
					_sdot_tmp.clear();

					s_cur = _s.back();
					sdot_cur = _sdot.back();
					beta_cur = determineAlphaBeta(s_cur, sdot_cur)(1);

					FI_SW = true;
					BI_SW = false;
				}
			}
		}

		// Step 3 : there exist two cases.
		//s_cur = _sf-1e-3;
		s_cur = _sf-0.0001;
		sdot_cur = _vf / _dqds(_sf-0.0001).norm();

		cout << "s_cur : " << s_cur << endl;
		cout << "sdot_cur : " << sdot_cur << endl;

		_s_tmp.push_front(s_cur);
		_sdot_tmp.push_front(sdot_cur);

		if (!swiPoint_swi) // case 1. when can't find switching point until final s
		{
			LOGIF(_sdot.back() > sdot_cur, "step 3 error(!swiPoint_swi) : s_end has to be smaller than _sdot.back().");
			
			while (s_cur >= _s.back())
			{
				s_BI_jk.push_back(s_cur);
				sd_BI_jk.push_back(sdot_cur);

				cout << "_s_back() : " << _s.back() << endl;
				alphabeta = determineAlphaBeta(s_cur, sdot_cur);
				alpha_cur = alphabeta(0);
				cout << "alpha_cur : " << alpha_cur << endl;
				cout << "beta_cur : " << beta_cur << endl;

				backwardIntegrate(s_cur, sdot_cur, alpha_cur);
				cout << "s_cur : " << s_cur << endl;
				cout << "sdot_cur : " << sdot_cur << endl;

				_s_tmp.push_front(s_cur);
				_sdot_tmp.push_front(sdot_cur);
			}

			_s_tmp.pop_front();
			_sdot_tmp.pop_front();

			sdot_cur = (sdot_cur - _sdot_tmp.front()) / (s_cur - _s_tmp.front()) * (_s.back() - s_cur) + sdot_cur;
			s_cur = _s.back();

			LOGIF(_sdot.back() > sdot_cur, "Backward intergration error:_sdot.back() has to be larger than sdot_cur.");
		}
		else // case 2. when forward integration reach final s
		{
			LOGIF(_sdot.back() > sdot_cur, "step 3 error(swiPoint_swi) : s_end has to be smaller than _sdot.back().");

			_s.pop_back();
			_sdot.pop_back();

			alphabeta = determineAlphaBeta(s_cur, sdot_cur);
			alpha_cur = alphabeta(0);

			Real delta = std::abs((s_cur-_s.back()));
			Real tmp = -2 * alpha_cur*delta + sdot_cur*sdot_cur;
			LOGIF((tmp > 0), "TOPP::backwardIntegrate error : the value has to be positive.");

			s_cur -= delta;
			sdot_cur = sqrt(tmp);
		}
		_s_tmp.push_front(s_cur);
		_sdot_tmp.push_front(sdot_cur);

		while (_sdot.back() > sdot_cur)
		{
			s_BI_jk.push_back(s_cur);
			sd_BI_jk.push_back(sdot_cur);

			_s.pop_back();
			_sdot.pop_back();
			alphabeta = determineAlphaBeta(s_cur, sdot_cur);
			alpha_cur = alphabeta(0);
			backwardIntegrate(s_cur, sdot_cur, alpha_cur);
			_s_tmp.push_front(s_cur);
			_sdot_tmp.push_front(sdot_cur);
		}

		// save backward information
		s_BI_jk.push_back(s_cur);
		sd_BI_jk.push_back(sdot_cur);
		//saveRealVector2txt(s_BI_jk, "C:/Users/crazy/Desktop/Time optimization/s_bsw.txt");
		//saveRealVector2txt(sd_BI_jk, "C:/Users/crazy/Desktop/Time optimization/sdot_bsw.txt");
		//saveRealVector2txt(s_BI_jk, "D:/jkkim/Documents/matlabTest/sBI.txt");
		//saveRealVector2txt(sd_BI_jk, "D:/jkkim/Documents/matlabTest/sdotBI.txt");

		_s_tmp.pop_front();
		_sdot_tmp.pop_front();

		_s.insert(_s.end(), _s_tmp.begin(), _s_tmp.end());
		_sdot.insert(_sdot.end(), _sdot_tmp.begin(), _sdot_tmp.end());
		_s_tmp.clear();
		_sdot_tmp.clear();

		// calculate tf and torque trajectory
		calculateFinalTime();
		//calculateTorqueTrajectory();

		std::list<Real>::iterator s_it;
		std::list<Real>::iterator sdot_it = _sdot.begin();

		for (s_it = _s.begin(); s_it != _s.end(); ++s_it)
		{
			s_jk.push_back(*(s_it));
			sdot_jk.push_back(*(sdot_it));
			sdot_it++;
		}
		//saveRealVector2txt(s_jk, "C:/Users/crazy/Desktop/Time optimization/s_result.txt");
		//saveRealVector2txt(sdot_jk, "C:/Users/crazy/Desktop/Time optimization/sdot_result.txt");
		//saveRealVector2txt(s_jk, "D:/jkkim/Documents/matlabTest/s_result.txt");
		//saveRealVector2txt(sdot_jk, "D:/jkkim/Documents/matlabTest/sdot_result.txt");
		saveRealVector2txt(s_jk, "C:/Users/ksh/Documents/MATLAB/s_result.txt");
		saveRealVector2txt(sdot_jk, "C:/Users/ksh/Documents/MATLAB/sdot_result.txt");

		cout << "trajectory generation finished." << endl;
	}

	const Real TOPP::getFinalTime() const
	{
		return _tf_result;
	}

	const MatrixX& TOPP::getTorqueTrajectory() const
	{
		return _torque_result;
	}

	////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////

	void TOPP::calcMVC()
	{
		Real ds = 1E-3;
		Real s = _si;
		Real sdot;
		while (true)
		{
			sdot = calculateMVCPoint(s);

			s_MVC_jk.push_back(s);
			sd_MVC_jk.push_back(sdot);

			//std::cout << s << std::endl;

			s += ds;
			if (s > _sf)
				break;
		}

	}
	void TOPP::calcSPs()
	{
		Real ds = 1E-3;
		Real s = _si;
		Real sd;
		int SPID;
		while (findNearestSwitchPoint(s))
		{
			std::cout << s << std::endl;

			s = _switchPoint[_switchPoint.size() - 1]._s;
			sd = _switchPoint[_switchPoint.size() - 1]._sdot;
			SPID = _switchPoint[_switchPoint.size() - 1]._id;

			s_SW_jk.push_back(s);
			sd_SW_jk.push_back(sd);
			SPID_SW_jk.push_back((Real)SPID);

			s += ds;
			if (s > _sf)
				break;
		}

	}
	void TOPP::saveRealVector2txt(std::vector<Real> in, std::string filename)
	{
		std::ofstream fout;
		fout.open(filename);

		for (unsigned int i = 0; i < in.size(); i++)
			fout << in[i] << std::endl;

		fout.close();

	}


	void TOPP::saveMVCandSP2txt()
	{
		//calcMVC();
		//saveRealVector2txt(s_MVC_jk, "C:/Users/crazy/Desktop/Time optimization/s.txt");
		//saveRealVector2txt(sd_MVC_jk, "C:/Users/crazy/Desktop/Time optimization/sdot.txt");


		////calcSPs();

		////saveRealVector2txt(s_SW_jk, "C:/Users/crazy/Desktop/Time optimization/s_sw.txt");
		////saveRealVector2txt(sd_SW_jk, "C:/Users/crazy/Desktop/Time optimization/sdot_sw.txt");

		//calcMVC();
		//saveRealVector2txt(s_MVC_jk, "D:/jkkim/Documents/matlabTest/sMVC.txt");
		//saveRealVector2txt(sd_MVC_jk, "D:/jkkim/Documents/matlabTest/sdotMVC.txt");

		//calcSPs();

		//saveRealVector2txt(s_SW_jk, "D:/jkkim/Documents/matlabTest/sSW.txt");
		//saveRealVector2txt(sd_SW_jk, "D:/jkkim/Documents/matlabTest/sdotSW.txt");


		calcSPs();
		saveRealVector2txt(s_SW_jk, "C:/Users/ksh/Documents/MATLAB/sSW.txt");
		saveRealVector2txt(sd_SW_jk, "C:/Users/ksh/Documents/MATLAB/sdotSW.txt");
		saveRealVector2txt(SPID_SW_jk, "C:/Users/ksh/Documents/MATLAB/spidSW.txt");

	}
}