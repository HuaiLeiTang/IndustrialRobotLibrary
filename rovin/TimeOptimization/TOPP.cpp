#include "TOPP.h"
#include <string>
#include <time.h>

using namespace std;
// a, b, c 미리 저장해서 쓰는거 생각하기 --> 계산 속도 향상 

namespace rovin {

	TOPP::TOPP(const MatrixX & q_data, const SerialOpenChainPtr & soc, const Real vi, const Real vf,
		const Real ds, const Real si, const Real sf, CONSTRAINT_TYPE constraintType)
	{
		// _q_data dimesion : dof * number of data
		_q_data = q_data;
		_soc = soc;
		_dof = _soc->getNumOfJoint();
		_state = _soc->makeState();
		_ds = ds;
		_vi = vi;
		_vf = vf;
		_si = si;
		_sf = sf;
		_constraintType = constraintType;

		// make b-spline
		makespline();

		// setting contraint
		settingconstraint();
	}

	TOPP::TOPP(const SerialOpenChainPtr & soc, const Real vi, const Real vf, const Real ds, const Real si, const Real sf, CONSTRAINT_TYPE constraintType)
	{
		_soc = soc;
		_dof = _soc->getNumOfJoint();
		_state = _soc->makeState();
		_ds = ds; _vi = vi; _vf = vf; _si = si; _sf = sf;
		_constraintType = constraintType;

		settingconstraint();
	}

	void TOPP::initialization()
	{
		_s.clear();
		_sdot.clear();
		_sddot.clear();
		_switchPoint.clear();
		_allMVCPoints.clear();
		_allMVCPointsFlag.clear();
	}

	void TOPP::setJointTrajectory(const MatrixX & q_data)
	{
		LOGIF((q_data.rows() == _dof), "TOPP::setJointTrajectory error : q_data size is wrong.");
		_q_data = q_data;
		makespline();
	}

	void TOPP::setSerialOpenChain(const SerialOpenChainPtr & soc)
	{
		_soc = soc;
		_dof = _soc->getNumOfJoint();
		LOGIF((_q_data.col(0).size() != _dof), "TOPP::setSerialOpenChain error : q_data size is wrong.");

		settingconstraint();
	}

	void TOPP::settingconstraint()
	{
		// insert torque constraint, velocity constraint, acceleration constraint
		_torqueConstraint.resize(_dof * 2);
		_velConstraint.resize(_dof * 2);
		_accConstraint.resize(_dof * 2);

		for (unsigned int i = 0; i < _dof; i++)
		{
			_torqueConstraint[i] = _soc->getMotorJointPtr(i)->getLimitTorqueUpper();
			_torqueConstraint[i + _dof] = _soc->getMotorJointPtr(i)->getLimitTorqueLower();

			_velConstraint[i] = _soc->getMotorJointPtr(i)->getLimitVelUpper();
			_velConstraint[i + _dof] = _soc->getMotorJointPtr(i)->getLimitVelLower();

			_accConstraint[i] = _soc->getMotorJointPtr(i)->getLimitAccUpper();
			_accConstraint[i + _dof] = _soc->getMotorJointPtr(i)->getLimitAccLower();
		}

		// calcuate number of constraints according to constraint type
		if (_constraintType == TORQUE)
		{
			_nconstraints = _torqueConstraint.size();
			_nconstraintsWithoutVel = _nconstraints;
		}
		else if (_constraintType == TORQUE_VEL)
		{
			_nconstraints = _torqueConstraint.size() + _velConstraint.size();
			_nconstraintsWithoutVel = _torqueConstraint.size();
		}
		else if (_constraintType == TORQUE_ACC)
		{
			_nconstraints = _torqueConstraint.size() + _accConstraint.size();
			_nconstraintsWithoutVel = _nconstraints;
		}
		else if (_constraintType == TORQUE_VEL_ACC)
		{
			_nconstraints = _torqueConstraint.size() + _velConstraint.size() + _accConstraint.size();
			_nconstraintsWithoutVel = _torqueConstraint.size() + _accConstraint.size();
		}
		else
			LOGIF(false, "TOPP::determineAlphaBeta error : contraint type is wrong.");
	}

	void TOPP::makespline()
	{
		unsigned int degreeOfBSpline = 3;
		unsigned int orderOfBSpline = degreeOfBSpline + 1;
		unsigned int MaxNumOfCP = 7;
		unsigned int numData = _q_data.cols();
		if (numData > MaxNumOfCP)
		{
			//cout << "Make b-spline using fitting" << endl;
			_q = BSplineFitting(_q_data, orderOfBSpline, MaxNumOfCP, _si, _sf);
		}
		else
		{
			//cout << "Make b-spline using interpolation" << endl;
			_q = BSplineInterpolation(_q_data, orderOfBSpline, _si, _sf);
		}
		_dqds = _q.derivative();
		_ddqdds = _dqds.derivative();
	}

	void TOPP::calculateAllMVCPoint()
	{
		int flag;
		Real s = _si;
		Real sdot = calculateMVCPoint(s, flag);
		Vector2 MVCPoint;

		while (s <= _sf)
		{
			MVCPoint(0) = s; MVCPoint(1) = sdot;
			_allMVCPoints.push_back(MVCPoint);
			_allMVCPointsFlag.push_back(flag);

			s += _ds;
			sdot = calculateMVCPoint(s, flag);
		}
		sdot = calculateMVCPoint(s - 1e-5, flag);
		MVCPoint(0) = s; MVCPoint(1) = sdot;
		_allMVCPoints.push_back(MVCPoint);
		_allMVCPointsFlag.push_back(flag);
	}

	void TOPP::calculateAllSwitchPoint()
	{
		Real s = _si;

		while (findNearestSwitchPoint(s))
		{
			s = _switchPoint[_switchPoint.size() - 1]._s;
			s += _ds;
			if (s > _sf)
				break;
		}
	}

	void TOPP::calculateA(Real s, VectorX& a)
	{
		VectorX q = _q(s), qs = _dqds(s);
		VectorX qdot = VectorX::Zero(q.size());
		//StatePtr state = _soc->makeState();
		_state->setJointStatePos(q); _state->setJointStateVel(qdot); _state->setJointStateAcc(qs);
		_soc->solveInverseDynamics(*_state);

		VectorX tau = _state->getJointStateTorque();

		qs.setZero();
		_state->setJointStatePos(q); _state->setJointStateVel(qdot); _state->setJointStateAcc(qs);
		_soc->solveInverseDynamics(*_state);

		VectorX tmp_a = tau - _state->getJointStateTorque();

		if (_constraintType == TORQUE || _constraintType == TORQUE_VEL)
		{
			a = VectorX(_torqueConstraint.size());
			for (unsigned int i = 0; i < _dof; i++)
			{
				a[i] = tmp_a[i];
				a[i + _dof] = -tmp_a[i];
			}
		}
		else if (_constraintType == TORQUE_ACC || _constraintType == TORQUE_VEL_ACC)
		{
			a = VectorX(_torqueConstraint.size() + _accConstraint.size());
			qs = _dqds(s);
			for (unsigned int i = 0; i < _dof; i++)
			{
				a[i] = tmp_a[i];
				a[i + _dof] = -tmp_a[i];
				a[i + _dof * 2] = qs(i);
				a[i + _dof * 3] = -qs(i);
			}
		}
		else
			LOGIF(false, "TOPP::calculateA error : contraint type is wrong.");
	}

	std::vector<VectorX> TOPP::calculateBandC(Real s)
	{
		VectorX q = _q(s), qdot = _dqds(s), qddot = _ddqdds(s);
		//StatePtr state = _soc->makeState();
		_state->setJointStatePos(q); _state->setJointStateVel(qdot); _state->setJointStateAcc(qddot);
		_soc->solveInverseDynamics(*_state);

		VectorX tau = _state->getJointStateTorque();

		_state->setJointStatePos(q);	_state->setJointStateVel(qdot.setZero()); _state->setJointStateAcc(qddot.setZero());
		_soc->solveInverseDynamics(*_state);

		VectorX tmp_c = _state->getJointStateTorque();
		VectorX tmp_b = tau - tmp_c;

		std::vector<VectorX> result;
		if (_constraintType == TORQUE || _constraintType == TORQUE_VEL)
		{
			VectorX b(_torqueConstraint.size()), c(_torqueConstraint.size());
			for (unsigned int i = 0; i < _dof; i++)
			{
				b[i] = tmp_b[i];
				b[i + _dof] = -tmp_b[i];
				c[i] = tmp_c[i];
				c[i + _dof] = -tmp_c[i];
			}
			for (unsigned int i = 0; i < _dof; i++)
			{
				c[i] -= _torqueConstraint[i];
				c[i + _dof] += _torqueConstraint[i + _dof];
			}
			result.push_back(b); result.push_back(c);
		}
		else if (_constraintType == TORQUE_ACC || _constraintType == TORQUE_VEL_ACC)
		{
			VectorX b(_torqueConstraint.size() + _accConstraint.size()), c(_torqueConstraint.size() + _accConstraint.size());
			VectorX qss = _ddqdds(s);
			for (unsigned int i = 0; i < _dof; i++)
			{
				b[i] = tmp_b[i];
				b[i + _dof] = -tmp_b[i];
				b[i + _dof * 2] = qss[i];
				b[i + _dof * 3] = -qss[i];
				c[i] = tmp_c[i];
				c[i + _dof] = -tmp_c[i];
				c[i + _dof * 2] = -_accConstraint[i];
				c[i + _dof * 3] = _accConstraint[i + _dof];
			}
			for (unsigned int i = 0; i < _dof; i++)
			{
				c[i] -= _torqueConstraint[i];
				c[i + _dof] += _torqueConstraint[i + _dof];
			}
			result.push_back(b); result.push_back(c);
		}
		else
			LOGIF(false, "TOPP::calculateBandC error : contraint type is wrong.");

		return result;
	}

	Vector2 TOPP::determineAlphaBeta(Real s, Real sdot, VectorX& a)
	{
		//if (s == _si)
		//	s += 1e-6;
		if (s == _sf)
			s -= 1e-6;

		//VectorX a;
		calculateA(s, a);

		VectorX q = _q(s), qdot = _dqds(s)*sdot, qddot = _ddqdds(s)*sdot*sdot;
		_state->setJointStatePos(q); _state->setJointStateVel(qdot); _state->setJointStateAcc(qddot);
		_soc->solveInverseDynamics(*_state);
		VectorX tau = _state->getJointStateTorque();

		VectorX left_vec;
		if (_constraintType == TORQUE || _constraintType == TORQUE_VEL)
		{
			left_vec = VectorX(_torqueConstraint.size());
			for (unsigned int i = 0; i < _dof; i++)
			{
				left_vec(i) = -tau(i) + _torqueConstraint(i);
				left_vec(i + _dof) = tau(i) - _torqueConstraint(i + _dof);
			}
		}
		else if (_constraintType == TORQUE_ACC || _constraintType == TORQUE_VEL_ACC)
		{
			left_vec = VectorX(_torqueConstraint.size() + _accConstraint.size());
			VectorX qss = _ddqdds(s)*sdot*sdot;
			for (unsigned int i = 0; i < _dof; i++)
			{
				left_vec(i) = -tau(i) + _torqueConstraint(i);
				left_vec(i + _dof) = tau(i) - _torqueConstraint(i + _dof);
				left_vec(i + _dof * 2) = -qss(i) + _accConstraint(i);
				left_vec(i + _dof * 3) = qss(i) - _accConstraint(i + _dof);
			}
		}
		else
			LOGIF(false, "TOPP::determineAlphaBeta error : contraint type is wrong.");

		Vector2 result; // result[0] is alpha, result[1] is beta
		result[0] = -std::numeric_limits<Real>::max();
		result[1] = std::numeric_limits<Real>::max();

		LOGIF(a.size() == left_vec.size(), "TOPP::determineAlphaBeta error : a & left_vec size are different.");

		for (int i = 0; i < left_vec.size(); i++)
		{
			Real tmp = left_vec(i) / a(i);
			if (a[i] > RealEps) // upper bound beta
			{
				if (tmp < result[1])
					result[1] = tmp;
			}
			else if (a[i] < -RealEps)// lower bound alpha
			{
				if (tmp > result[0])
					result[0] = tmp;
			}
		}
		return result;
	}

	void TOPP::determineVelminmax(Real s, Vector2& minmax)
	{
		VectorX qs = _dqds(s), qs_vec(_dof * 2), left_vec(_dof * 2);
		for (unsigned int i = 0; i < _dof; i++)
		{
			qs_vec(i) = qs(i);
			qs_vec(i + _dof) = -qs(i);
			left_vec(i) = _velConstraint(i);
			left_vec(i + _dof) = -_velConstraint(i + _dof);
		}
		minmax(0) = -std::numeric_limits<Real>::max();
		minmax(1) = std::numeric_limits<Real>::max();
		for (unsigned int i = 0; i < _dof * 2; i++)
		{
			Real tmp = left_vec(i) / qs_vec(i);
			if (qs_vec(i) > RealEps) // upper bound beta
			{
				if (tmp < minmax(1))
					minmax(1) = tmp;
			}
			else if (qs_vec(i) > -RealEps)// lower bound alpha
			{
				if (tmp > minmax(0))
					minmax(0) = tmp;
			}
		}
		minmax(0) = std::max(minmax(0), 0.0);
	}

	Real TOPP::calculateMVCPoint(Real s, int& flag)
	{
		Real sdot_MVC = std::numeric_limits<Real>::max();

		VectorX a;
		Vector2 alphabeta = determineAlphaBeta(s, 0, a);
		if (alphabeta[0] > alphabeta[1])
			return 0;

		//calculateA(s, a);

		std::vector<VectorX> BandC = calculateBandC(s);
		VectorX b = BandC[0], c = BandC[1];

		LOGIF(((a.size() == b.size()) && (a.size() == c.size()) && (b.size() == c.size())), "TOPP::calulateMVCPoint error : a, b, c vector size is wrong.");

		Real sdot_VelC = std::numeric_limits<Real>::max();
		Vector2 minmax;
		if (_constraintType == TORQUE_VEL || _constraintType == TORQUE_VEL_ACC)
		{
			determineVelminmax(s, minmax);
			sdot_VelC = minmax(1);
		}

		unsigned int kk;
		unsigned int mm;

		for (unsigned int k = 0; k < _nconstraintsWithoutVel; k++)
		{
			for (unsigned int m = k + 1; m < _nconstraintsWithoutVel; m++)
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
		if (sdot_MVC <= sdot_VelC)
		{
			flag = 2;
			return sdot_MVC;
		}
		else
		{
			flag = 1;
			return sdot_VelC;
		}
	}

	Real TOPP::calculateMVCPointExclude(Real s, int iExclude, int& flag)
	{
		Real sdot_MVC = std::numeric_limits<Real>::max();

		VectorX a;
		Vector2 alphabeta = determineAlphaBeta(s, 0, a);
		if (alphabeta[0] > alphabeta[1])
			return 0;

		//calculateA(s, a);
		std::vector<VectorX> BandC = calculateBandC(s);
		VectorX b = BandC[0], c = BandC[1];

		LOGIF(((a.size() == b.size()) && (a.size() == c.size()) && (b.size() == c.size())), "TOPP::calulateMVCPoint error : a, b, c vector size is wrong.");

		Real sdot_VelC = std::numeric_limits<Real>::max();
		Vector2 minmax;
		if (_constraintType == TORQUE_VEL || _constraintType == TORQUE_VEL_ACC)
		{
			determineVelminmax(s, minmax);
			sdot_VelC = minmax(1);
		}

		for (unsigned int k = 0; k < _nconstraintsWithoutVel; k++)
		{
			for (unsigned int m = k + 1; m < _nconstraintsWithoutVel; m++)
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
		if (sdot_MVC <= sdot_VelC)
		{
			flag = 2;
			return sdot_MVC;
		}
		else
		{
			flag = 1;
			return sdot_VelC;
		}
	}

	void TOPP::calculateFinalTime()
	{
		_t = VectorX(_s.size());
		std::list<Real>::iterator it_sdot = ++(_sdot.begin());
		Real s_k = _s.front(), sdot_k = _sdot.front();
		Real s_k1, sdot_k1;
		Real sum = 0;
		unsigned int cnt = 0;
		_t(cnt) = sum; cnt++;
		for (std::list<Real>::iterator it_s = ++(_s.begin()); it_s != (_s.end()); ++it_s)
		{
			s_k1 = *(it_s); sdot_k1 = *(it_sdot);
			sum += 2 * (s_k1 - s_k) / (sdot_k1 + sdot_k);
			_t(cnt) = sum;
			s_k = s_k1; sdot_k = sdot_k1;
			it_sdot++; cnt++;
		}
		_tf_result = sum;

		//std::cout << _t << std::endl;
	}

	void TOPP::calculateTorqueTrajectory()
	{
		LOGIF((_t.size() != 0), "calculateTorqueTrajectory error : you should calculate final time.");

		_torque_result.resize(_dof, _s.size() - 1);
		std::list<Real>::iterator it_sdot = _sdot.begin();
		std::list<Real>::iterator it_sddot = _sddot.begin();
		unsigned int cnt = 0;
		VectorX q, qdot, qddot;
		VectorX qs;

		for (std::list<Real>::iterator it_s = _s.begin(); it_s != --(_s.end()); it_s++)
		{
			q = _q(*it_s); qs = _dqds(*it_s);
			qdot = (*it_sdot) * qs;
			qddot = (*it_sddot) * qs + (*it_sdot) * (*it_sdot) * _ddqdds(*it_s);

			_state->setJointStatePos(q);
			_state->setJointStateVel(qdot);
			_state->setJointStateAcc(qddot);
			_soc->solveInverseDynamics(*_state);

			_torque_result.col(cnt) = _state->getJointStateTorque();

			it_sdot++; it_sddot++;
			cnt++;
		}
	}


	void TOPP::forwardIntegrate(Real & s, Real & sdot, Real sddot)
	{
		Real tmp = 2 * sddot*_ds + sdot*sdot;
		LOGIF((tmp > 0), "TOPP::forwardIntegrate error : the value has to be positive.");
		s = s + _ds;
		sdot = sqrt(tmp);
	}

	void TOPP::backwardIntegrate(Real & s, Real & sdot, Real sddot)
	{
		Real tmp = -2 * sddot*_ds + sdot*sdot;
		LOGIF((tmp > 0), "TOPP::backwardIntegrate error : the value has to be positive.");
		s = s - _ds;
		sdot = sqrt(tmp);
	}

	bool TOPP::findNearestSwitchPoint(Real s)
	{
		int flag, flag_cur, flag_next;


		//Real ds = 0.0005;
		Real ds = _ds;
		Real s_bef = s;
		Real sdot_bef = calculateMVCPoint(s_bef, flag);
		Real s_cur = s_bef + ds;
		Real sdot_cur = calculateMVCPoint(s_cur, flag_cur);
		Real s_next = s_bef + 2 * ds;
		Real sdot_next = calculateMVCPoint(s_next, flag_next);

		// end criterion
		if (s_next >= _sf)
			false;

		// 원래 topp 코드에는 tan_bef/tan_cur 로 discontinous 포인트 추가하는데 그거 이해 안감
		// 안넣어도 될것같은데 그럴거면 tan_bef/tan_cur 따로 구할 필요 없음
		//Real tan_bef = (sdot_cur - sdot_bef) / ds;
		//Real tan_cur = (sdot_next - sdot_cur) / ds;
		VectorX a_bef;
		VectorX a_cur;
		Real diff_bef = determineAlphaBeta(s_bef, sdot_bef, a_bef)[0] / sdot_bef - (sdot_cur - sdot_bef) / ds;
		Real diff_cur = determineAlphaBeta(s_cur, sdot_cur, a_cur)[0] / sdot_cur - (sdot_next - sdot_cur) / ds;

		//calculateA(s_bef, a_bef);
		std::vector<VectorX> bc_bef = calculateBandC(s_bef);
		//calculateA(s_cur, a_cur);
		std::vector<VectorX> bc_cur = calculateBandC(s_cur);


		bool singularFound = false;
		Real s_sing, a_sing, b_sing, c_sing, sdot_star, sdot_plus, s_save; ///> variables for finding singular points
		Real sdot_min = std::numeric_limits<Real>::max();
		int idx_save;

		while (true)
		{

			if (flag_cur == 2 && flag_next == 2)
			{

				// step 1: singular point check
				// include calculating \lambda
				for (unsigned int i = 0; i < _nconstraintsWithoutVel; i++) // number of inequality
				{
					if (a_bef[i] * a_cur[i] <= 0)
					{
						// 딱 0이 아니면 원래 코드처럼 선형보간 해서 더 좋은 s를 찾는 과정이 필요할듯..
						// 여기 말고 위아래 step 2, 3에도 비슷한 알고리즘 필요할듯
						Real adiff = a_cur[i] - a_bef[i];
						if (std::abs(adiff) > 1E-10)
						{
							// calculate '_sing' variables as internal dividing points
							// abc도 선형보간으로?? 새로 구한 s_sing 에서 구하는게 아니고??
							Real tmpalpha = -a_bef[i] / adiff;
							s_sing = tmpalpha*s_cur + (1 - tmpalpha)*s_bef;
							VectorX tmp_sing;
							calculateA(s_sing, tmp_sing);
							a_sing = tmp_sing[i];
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
							sdot_plus = calculateMVCPointExclude(s_sing, i, flag);

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
			}
			if (singularFound)
			{
				Real diffeps = 1E-5;
				Real ap, bp, cp; ///> p for prime, differentiated to s
				VectorX tmp1, tmp2;
				calculateA(s_save + diffeps, tmp1); calculateA(s_save, tmp2);
				ap = (tmp1(idx_save) - tmp2(idx_save)) / diffeps;
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
			flag_cur = flag_next;
			s_next += ds;
			sdot_next = calculateMVCPoint(s_next, flag_next);

			a_bef = a_cur;

			diff_bef = diff_cur;
			diff_cur = determineAlphaBeta(s_cur, sdot_cur, a_cur)[0] / sdot_cur - (sdot_next - sdot_cur) / ds;

			bc_bef = bc_cur;
			//calculateA(s_cur, a_cur);
			bc_cur = calculateBandC(s_cur);


			// end criterion
			if (s_next >= _sf)
				return false;
		}
	}

	unsigned int TOPP::forward(Real& s_cur, Real& sdot_cur)
	{
		Real beta_cur, sdot_cur_MVC;
		unsigned int idx;
		int flag;

		while (std::abs(s_cur - _sf) > _ds*0.1)
		{
			beta_cur = determineAlphaBeta(s_cur, sdot_cur)(1);
			forwardIntegrate(s_cur, sdot_cur, beta_cur);
			if (sdot_cur < 1e-5)
				return 0;
			idx = (unsigned int)round(s_cur / _ds);
			sdot_cur_MVC = _allMVCPoints[idx](1);
			flag = _allMVCPointsFlag[idx];

			if (sdot_cur > sdot_cur_MVC)
			{
				sdot_cur = sdot_cur_MVC;
				_s.push_back(s_cur); _sdot.push_back(sdot_cur);
				_sddot.push_back(beta_cur);

				if (flag == 2)
					return 1;
				if (flag == 1)
					return 2;
			}
			_s.push_back(s_cur); _sdot.push_back(sdot_cur);
			_sddot.push_back(beta_cur);
		}

		_s.pop_back(); _sdot.pop_back();
		s_cur = _sf; sdot_cur = _vf / _dqds(_sf - 1e-4).norm();
		_s.push_back(s_cur); _sdot.push_back(sdot_cur);

		return 3;
	}

	unsigned int TOPP::forwardVel(Real& s_cur, Real& sdot_cur)
	{
		Vector2 alphabeta;
		Real s_next, sdot_next, alpha_cur, beta_cur, slope;
		unsigned int idx;
		int flag;

		while (true)
		{
			s_next = s_cur + _ds;

			idx = (unsigned int)round(s_next / _ds);
			sdot_next = _allMVCPoints[idx](1);
			flag = _allMVCPointsFlag[idx];

			alphabeta = determineAlphaBeta(s_cur, sdot_cur);
			alpha_cur = alphabeta(0); beta_cur = alphabeta(1);
			slope = (sdot_next - sdot_cur) / (s_next - s_cur);

			if (flag == 2)
			{
				//_sddot.pop_back();
				s_cur = s_next; sdot_cur = sdot_next;
				_s.push_back(s_cur); _sdot.push_back(sdot_cur); _sddot.push_back(slope*sdot_cur);
				return 2;
			}

			if (std::abs(s_next - _sf) < _ds*0.1)
			{
				s_cur = _sf; sdot_cur = _vf / _dqds(_sf - 1e-4).norm();
				_s.push_back(s_cur); _sdot.push_back(sdot_cur); _sddot.push_back(slope*sdot_cur);
				return 0;
			}
					
			if (beta_cur >= slope && slope >= alpha_cur)
			{
				s_cur = s_next; sdot_cur = sdot_next;
				_s.push_back(s_cur); _sdot.push_back(sdot_cur);
				_sddot.push_back(slope*sdot_cur);
			}
			else if (beta_cur < slope)
				return 1;
			else
				break;
		}

		while (true)
		{	
			s_next = s_cur + _ds;

			idx = (unsigned int)round(s_next / _ds);
			sdot_next = _allMVCPoints[idx](1);
			flag = _allMVCPointsFlag[idx];

			alphabeta = determineAlphaBeta(s_cur, sdot_cur);
			alpha_cur = alphabeta(0); beta_cur = alphabeta(1);
			slope = (sdot_next - sdot_cur) / (s_next - s_cur);

			if (flag == 2)
			{
				s_cur = s_next; sdot_cur = sdot_next;
				return 2;
			}

			if (std::abs(s_next - _sf) < _ds*0.1)
			{
				s_cur = _sf; sdot_cur = _vf / _dqds(_sf - 1e-4).norm();
				return 0;
			}
			
			if (beta_cur >= slope && slope >= alpha_cur)
			{
				s_cur = s_next; sdot_cur = sdot_next;
				return 3;
			}
			s_cur = s_next; sdot_cur = sdot_next;
		}
	}

	void TOPP::backward(Real& s_cur, Real& sdot_cur)
	{
		Real alpha_cur;
		Real s_back = _s.back();
		std::list<Real> _s_tmp;
		std::list<Real> _sdot_tmp;
		std::list<Real> _sddot_tmp;

		while (true)
		{
			if (std::abs(s_cur - s_back) <= _ds*0.1)
			{
				//LOGIF(sdot_cur < _sdot.back(), "backward error : sdot_cur must be smaller than _sdot.back().");
				alpha_cur = determineAlphaBeta(s_cur, sdot_cur)(0);
				_s_tmp.push_front(s_cur); _sdot_tmp.push_front(sdot_cur); _sddot_tmp.push_front(alpha_cur);
				_s.pop_back(); _sdot.pop_back(); 
				//_sddot.pop_back();

				while (true)
				{					
					backwardIntegrate(s_cur, sdot_cur, alpha_cur);
					alpha_cur = determineAlphaBeta(s_cur, sdot_cur)(0);
					if (sdot_cur > _sdot.back())
						break;

					_s_tmp.push_front(s_cur); _sdot_tmp.push_front(sdot_cur); _sddot_tmp.push_front(alpha_cur);
					_s.pop_back(); _sdot.pop_back(); _sddot.pop_back();
				}
				break;
			}
			_s_tmp.push_front(s_cur); _sdot_tmp.push_front(sdot_cur);

			alpha_cur = determineAlphaBeta(s_cur, sdot_cur)(0);
			_sddot_tmp.push_front(alpha_cur);
			backwardIntegrate(s_cur, sdot_cur, alpha_cur);
		}
		_sddot_tmp.pop_back();
		s_cur = _s_tmp.back(); sdot_cur = _sdot_tmp.back();
		_s.insert(_s.end(), _s_tmp.begin(), _s_tmp.end());
		_sdot.insert(_sdot.end(), _sdot_tmp.begin(), _sdot_tmp.end());
		_sddot.insert(_sddot.end(), _sddot_tmp.begin(), _sddot_tmp.end());
	}

	void TOPP::switchpointfunction(Real& s_cur, Real& sdot_cur, bool& singularPoint_swi)
	{
		bool swiPoint_swi;
		swiPoint_swi = findNearestSwitchPoint(s_cur);
		if (!swiPoint_swi)
		{
			s_cur = _sf; 
			sdot_cur = _vf / _dqds(_sf - 1e-4).norm();
		}
		else
		{
			s_cur = _switchPoint[_switchPoint.size() - 1]._s;
			sdot_cur = _switchPoint[_switchPoint.size() - 1]._sdot;
			s_cur = round(s_cur / _ds) * _ds;										
			if (_switchPoint[_switchPoint.size() - 1]._id == SwitchPoint::SINGULAR)
				singularPoint_swi = true;
		}
	}

	bool TOPP::generateTrajectory()
	{
		initialization();
		calculateAllMVCPoint();

		Real s_cur = 0;
		Real sdot_cur = _vi / _dqds(1e-5).norm();
		Vector2 alphabeta = determineAlphaBeta(s_cur, sdot_cur);
		Real alpha_cur = alphabeta(0);
		Real beta_cur = alphabeta(1);
		unsigned int numOfSPInt = 3;
		unsigned int swi_forward;
		unsigned int swi_forwardVel = 1;
		bool singularPoint_swi = false;

		_s.push_back(s_cur); _sdot.push_back(sdot_cur);

		while (true)
		{
			if (swi_forwardVel == 1)
			{
				if (singularPoint_swi)
				{
					Real lambda = _switchPoint[_switchPoint.size() - 1]._lambda;
					for (unsigned int i = 0; i < numOfSPInt*2; i++)
					{
						s_cur += _ds;
						sdot_cur += lambda * _ds;
						_s.push_back(s_cur); _sdot.push_back(sdot_cur);
						_sddot.push_back(lambda*sdot_cur);
					}
					singularPoint_swi = false;
				}
				swi_forward = forward(s_cur, sdot_cur);
			}

			if (swi_forward == 0)
			{
				LOG("forward fail, TOPP fail.");
				return false;
			}
			else if (swi_forward == 1)
			{
				swi_forwardVel = 0; 
				switchpointfunction(s_cur, sdot_cur, singularPoint_swi);
			}
			else if (swi_forward == 2)
				swi_forwardVel = forwardVel(s_cur, sdot_cur);


			if (swi_forwardVel == 2)
			{
				swi_forwardVel = 0;
				switchpointfunction(s_cur, sdot_cur, singularPoint_swi);
			}
			if (swi_forwardVel == 0 || swi_forwardVel == 3 || swi_forward == 3)
			{
				if (singularPoint_swi)
				{
					Real lambda = _switchPoint[_switchPoint.size() - 1]._lambda;
					for (unsigned int i = 0; i < numOfSPInt; i++)
					{
						s_cur -= _ds;
						sdot_cur -= lambda * _ds;
					}
				}
				backward(s_cur, sdot_cur);

				if (swi_forwardVel == 0)
					swi_forwardVel = 1;
				else
				{
					swi_forwardVel = 3;
					swi_forward = 2;
				}
			}
			if (std::abs(_s.back() - _sf) < _ds*0.1)
				break;
		}

		std::cout << "trajectory generation finished." << std::endl;
		return true;
	}
}