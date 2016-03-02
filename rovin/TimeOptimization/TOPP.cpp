#include "TOPP.h"
#include <string>

using namespace std;

namespace rovin {

	TOPP::TOPP(const MatrixX & q_data, const SerialOpenChainPtr & soc, 
		const Real ds, const Real vi, const Real vf, const Real si, const Real sf, CONSTRAINT_TYPE constraintType)
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
		_constraintType = constraintType;

		// insert torque constraint, velocity constraint, acceleration constraint
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
			_accConstraint[i + _dof] = _soc->getMotorJointPtr(i)->getLimitAccLower();
		}

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
		{
			LOGIF(false, "TOPP::determineAlphaBeta error : contraint type is wrong.");
		}
	}

	void TOPP::initialization()
	{

	}

	void TOPP::calculateA(Real s, VectorX& a)
	{
		VectorX q = _q(s), qs = _dqds(s);
		VectorX qdot = VectorX::Zero(q.size());
		StatePtr state = _soc->makeState();
		state->setJointStatePos(q);	state->setJointStateVel(qdot); state->setJointStateAcc(qs);
		_soc->solveInverseDynamics(*state);
		
		VectorX tau = state->getJointStateTorque();

		qs.setZero();
		state->setJointStatePos(q);	state->setJointStateVel(qdot); state->setJointStateAcc(qs);
		_soc->solveInverseDynamics(*state);

		VectorX tmp_a = tau - state->getJointStateTorque();

		if (_constraintType == TORQUE || _constraintType == TORQUE_VEL)
		{
			a = VectorX(_torqueConstraint.size());
			for (int i = 0; i < _dof; i++)
			{
				a[i] = tmp_a[i];
				a[i + _dof] = -tmp_a[i];
			}
		}
		else if (_constraintType == TORQUE_ACC || _constraintType == TORQUE_VEL_ACC)
		{
			a = VectorX(_torqueConstraint.size() + _accConstraint.size());
			qs = _dqds(s);
			for (int i = 0; i < _dof; i++)
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
		StatePtr state = _soc->makeState();
		state->setJointStatePos(q); state->setJointStateVel(qdot); state->setJointStateAcc(qddot);
		_soc->solveInverseDynamics(*state);

		VectorX tau = state->getJointStateTorque();

		state->setJointStatePos(q);	state->setJointStateVel(qdot.setZero()); state->setJointStateAcc(qddot.setZero());
		_soc->solveInverseDynamics(*state);

		VectorX tmp_c = state->getJointStateTorque();
		VectorX tmp_b = tau - tmp_c;
		
		std::vector<VectorX> result;
		if (_constraintType == TORQUE || _constraintType == TORQUE_VEL)
		{
			VectorX b(_torqueConstraint.size()), c(_torqueConstraint.size());
			for (int i = 0; i < _dof; i++)
			{
				b[i] = tmp_b[i];
				b[i + _dof] = -tmp_b[i];
				c[i] = tmp_c[i];
				c[i + _dof] = -tmp_c[i];
			}
			for (int i = 0; i < _dof; i++)
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
			for (int i = 0; i < _dof; i++)
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
			for (int i = 0; i < _dof; i++)
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

	Vector2 TOPP::determineAlphaBeta(Real s, Real sdot)
	{

		VectorX a;
		calculateA(s, a);
			
		VectorX q = _q(s), qdot = _dqds(s)*sdot, qddot = _ddqdds(s)*sdot*sdot;
		StatePtr state = _soc->makeState();
		state->setJointStatePos(q); state->setJointStateVel(qdot); state->setJointStateAcc(qddot);
		_soc->solveInverseDynamics(*state);
		VectorX tau = state->getJointStateTorque();

		VectorX left_vec;
		if (_constraintType == TORQUE || _constraintType == TORQUE_VEL)
		{
			left_vec = VectorX(_torqueConstraint.size());
			for (int i = 0; i < _dof; i++)
			{
				left_vec(i) = -tau(i) + _torqueConstraint(i);
				left_vec(i + _dof) = tau(i) - _torqueConstraint(i + _dof);
			}
		}
		else if (_constraintType == TORQUE_ACC || _constraintType == TORQUE_VEL_ACC)
		{
			left_vec = VectorX(_torqueConstraint.size()+ _accConstraint.size());
			VectorX qss = _ddqdds(s)*sdot*sdot;
			for (int i = 0; i < _dof; i++)
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
		for (int i = 0; i < _dof; i++)
		{
			qs_vec(i) = qs(i);
			qs_vec(i + _dof) = -qs(i);
			left_vec(i) = _velConstraint(i);
			left_vec(i + _dof) = -_velConstraint(i + _dof);
		}
		minmax(0) = -std::numeric_limits<Real>::max();
		minmax(1) = std::numeric_limits<Real>::max();
		for (int i = 0; i < _dof * 2; i++)
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

		Vector2 alphabeta = determineAlphaBeta(s, 0);
		if (alphabeta[0] > alphabeta[1])
			return 0;

		VectorX a;
		calculateA(s, a);
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

		for (int k = 0; k < _nconstraintsWithoutVel; k++)
		{
			for (int m = k + 1; m < _nconstraintsWithoutVel; m++)
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

		Vector2 alphabeta = determineAlphaBeta(s, 0);
		if (alphabeta[0] > alphabeta[1])
			return 0;

		VectorX a;
		calculateA(s, a);
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

		for (int k = 0; k < _nconstraintsWithoutVel; k++)
		{
			for (int m = k + 1; m < _nconstraintsWithoutVel; m++)
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
		int flag;

		Real ds = 0.0005;
		Real s_bef = s;
		Real sdot_bef = calculateMVCPoint(s_bef, flag);
		Real s_cur = s_bef + ds;
		Real sdot_cur = calculateMVCPoint(s_cur, flag);
		Real s_next = s_bef + 2 * ds;
		Real sdot_next = calculateMVCPoint(s_next, flag);

		// end criterion
		if (s_next >= _sf)
			false;

		// 원래 topp 코드에는 tan_bef/tan_cur 로 discontinous 포인트 추가하는데 그거 이해 안감
		// 안넣어도 될것같은데 그럴거면 tan_bef/tan_cur 따로 구할 필요 없음
		//Real tan_bef = (sdot_cur - sdot_bef) / ds;
		//Real tan_cur = (sdot_next - sdot_cur) / ds;
		Real diff_bef = determineAlphaBeta(s_bef, sdot_bef)[0] / sdot_bef - (sdot_cur - sdot_bef) / ds;
		Real diff_cur = determineAlphaBeta(s_cur, sdot_cur)[0] / sdot_cur - (sdot_next - sdot_cur) / ds;

		VectorX a_bef;
		calculateA(s_bef, a_bef);
		std::vector<VectorX> bc_bef = calculateBandC(s_bef);
		VectorX a_cur;
		calculateA(s_cur, a_cur);
		std::vector<VectorX> bc_cur = calculateBandC(s_cur);


		bool singularFound = false;
		Real s_sing, a_sing, b_sing, c_sing, sdot_star, sdot_plus, s_save; ///> variables for finding singular points
		Real sdot_min = std::numeric_limits<Real>::max();
		int idx_save;

		while (true)
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
			s_next += ds;
			sdot_next = calculateMVCPoint(s_next, flag);

			diff_bef = diff_cur;
			diff_cur = determineAlphaBeta(s_cur, sdot_cur)[0] / sdot_cur - (sdot_next - sdot_cur) / ds;

			a_bef = a_cur;
			bc_bef = bc_cur;
			calculateA(s_cur, a_cur);
			bc_cur = calculateBandC(s_cur);

			// end criterion
			if (s_next >= _sf)
				return false;
		}
	}

	void TOPP::generateTrajectory()
	{
		Real s_cur = 0, sdot_cur = _vi / _dqds(0.0001).norm(), sdot_MVC;
		_s.push_back(s_cur); _sdot.push_back(sdot_cur);
		
		Vector2 alphabeta = determineAlphaBeta(s_cur, sdot_cur);
		Real alpha_cur = alphabeta(0);
		Real beta_cur = alphabeta(1);

		std::list<Real> _s_tmp;
		std::list<Real> _sdot_tmp;

		int flag; ///< 1 : velocity contraint, 2 : acc-torque contraint
		bool FI_SW = true, BI_SW = false, I_SW = true;
		bool swiPoint_swi;
		unsigned int numOfSPInt = 3; ///< singular point integration number
		

		// Step 1 & 2 : Forward & backward integration
		while (I_SW)
		{
			while (FI_SW) ///< Forward integration
			{
				forwardIntegrate(s_cur, sdot_cur, beta_cur);
				_s.push_back(s_cur); _sdot.push_back(sdot_cur);
				
				beta_cur = determineAlphaBeta(s_cur, sdot_cur)(1);
				sdot_MVC = calculateMVCPoint(s_cur, flag);

				if(sdot_cur >= sdot_MVC)
				{
					sdot_cur = sdot_MVC; _sdot.pop_back(); _sdot.push_back(sdot_cur);

					if (flag == 1) ///< velocity contraint case
					{
						Real s_tmp_beta, sdot_tmp_beta;
						Real s_tmp_alpha, sdot_tmp_alpha;
						Real s_next, sdot_next;
						while (true)
						{
							s_tmp_beta = s_cur; sdot_tmp_beta = sdot_cur;
							s_tmp_alpha = s_cur; sdot_tmp_alpha = sdot_cur;
							s_next = s_cur + _ds;
							sdot_next = calculateMVCPoint(s_next, flag);
							
							if (flag == 2)
							{
								swiPoint_swi = findNearestSwitchPoint(s_cur);
								if (!swiPoint_swi)
								{
									FI_SW = false; BI_SW = false; I_SW = false;
								}
								else
								{
									s_cur = _switchPoint[_switchPoint.size() - 1]._s;
									sdot_cur = _switchPoint[_switchPoint.size() - 1]._sdot;
									_s_tmp.push_front(s_cur); _sdot_tmp.push_front(sdot_cur);
									if (_switchPoint[_switchPoint.size() - 1]._id == SwitchPoint::SINGULAR)
									{
										Real lambda = _switchPoint[_switchPoint.size() - 1]._lambda;
										for (int i = 0; i < numOfSPInt; i++)
										{
											s_cur -= _ds;
											sdot_cur -= lambda * _ds;
											_s_tmp.push_front(s_cur);
											_sdot_tmp.push_front(sdot_cur);
										}
									}
									FI_SW = false; BI_SW = true;
								}
								break;
							}

							alphabeta = determineAlphaBeta(s_cur, sdot_cur);
							alpha_cur = alphabeta(0);
							beta_cur = alphabeta(1);
							forwardIntegrate(s_tmp_beta, sdot_tmp_beta, beta_cur);
							forwardIntegrate(s_tmp_alpha, sdot_tmp_alpha, alpha_cur);
							if (sdot_tmp_beta >= sdot_next && sdot_next >= sdot_tmp_alpha)
							{
								s_cur = s_next; sdot_cur = sdot_next;
								_s.push_back(s_cur); _sdot.push_back(sdot_cur);
							}
							else if (sdot_tmp_beta < sdot_next)
							{
								FI_SW = true; BI_SW = false;
								break;
							}
							else if (sdot_tmp_alpha > sdot_next)
							{
								while (sdot_tmp_alpha > sdot_next)
								{
									s_tmp_alpha = s_cur; sdot_tmp_alpha = sdot_cur;
									s_next = s_cur + _ds;
									sdot_next = calculateMVCPoint(s_next, flag);

									if (flag == 2)
									{
										swiPoint_swi = findNearestSwitchPoint(s_cur);
										if (!swiPoint_swi)
										{
											FI_SW = false; BI_SW = false; I_SW = false;
										}
										else
										{
											s_cur = _switchPoint[_switchPoint.size() - 1]._s;
											sdot_cur = _switchPoint[_switchPoint.size() - 1]._sdot;
											_s_tmp.push_front(s_cur); _sdot_tmp.push_front(sdot_cur);
											if (_switchPoint[_switchPoint.size() - 1]._id == SwitchPoint::SINGULAR)
											{
												Real lambda = _switchPoint[_switchPoint.size() - 1]._lambda;
												for (int i = 0; i < numOfSPInt; i++)
												{
													s_cur -= _ds;
													sdot_cur -= lambda * _ds;
													_s_tmp.push_front(s_cur);
													_sdot_tmp.push_front(sdot_cur);
												}
											}
											FI_SW = false; BI_SW = true;
										}
										break;
									}

									alpha_cur = determineAlphaBeta(s_cur, sdot_cur)(0);
									forwardIntegrate(s_tmp_alpha, sdot_tmp_alpha, alpha_cur);
									s_cur = s_next; sdot_cur = sdot_next;
								}
								s_cur -= _ds;
								sdot_cur = calculateMVCPoint(s_cur, flag);

								FI_SW = false; BI_SW = true;
								break;
							}
							else if (s_cur > _sf)
							{
								FI_SW = true; BI_SW = false; I_SW = false;
								break;
							}
						}
					}
					else if (flag == 2) /// acc-torque constraint case
					{
						swiPoint_swi = findNearestSwitchPoint(s_cur);
						if (!swiPoint_swi) /// if swtich point doesn't exist until s_end --> go to step 3
						{
							FI_SW = false; BI_SW = false; I_SW = false;
						}
						else /// if switch point exist --> go to backward integration
						{
							s_cur = _switchPoint[_switchPoint.size() - 1]._s;
							sdot_cur = _switchPoint[_switchPoint.size() - 1]._sdot;
							_s_tmp.push_front(s_cur); _sdot_tmp.push_front(sdot_cur);
							if (_switchPoint[_switchPoint.size() - 1]._id == SwitchPoint::SINGULAR)
							{
								Real lambda = _switchPoint[_switchPoint.size() - 1]._lambda;
								for (int i = 0; i < numOfSPInt; i++)
								{
									s_cur -= _ds;
									sdot_cur -= lambda * _ds;
									_s_tmp.push_front(s_cur);
									_sdot_tmp.push_front(sdot_cur);
								}
							}
							FI_SW = false; BI_SW = true;
						}
					}
				}
				else if (s_cur > _sf) ///< go to step 3, case (a), s_cur 가 무조건 s_end 넘어갔을 때!!
				{
					FI_SW = false; BI_SW = false; I_SW = false;
				}
				else if (sdot_cur < 1e-10) ///< Debugging 요소, case (a)
				{
					FI_SW = false; BI_SW = false; I_SW = false;
				}
			}

			while (BI_SW) ///< Backward integration
			{

				alphabeta = determineAlphaBeta(s_cur, sdot_cur);
				alpha_cur = alphabeta(0);
				beta_cur = alphabeta(1);

				backwardIntegrate(s_cur, sdot_cur, alpha_cur);
				_s_tmp.push_front(s_cur);
				_sdot_tmp.push_front(sdot_cur);

				if (s_cur <= _s.back())
				{
					_s_tmp.pop_front();	_sdot_tmp.pop_front();
					sdot_cur = (sdot_cur - _sdot_tmp.front()) / (s_cur - _s_tmp.front()) * (_s.back() - s_cur) + sdot_cur;
					s_cur = _s.back();

					LOGIF(_sdot.back() > sdot_cur,"Backward intergration error:_sdot.back() has to be larger than sdot_cur.");
					_s_tmp.push_front(s_cur); _sdot_tmp.push_front(sdot_cur);

					while (_sdot.back() > sdot_cur)
					{
						_s.pop_back();
						_sdot.pop_back();
						alphabeta = determineAlphaBeta(s_cur, sdot_cur);
						alpha_cur = alphabeta(0);

						backwardIntegrate(s_cur, sdot_cur, alpha_cur); ///< update s_cur, sdot_cur
						_s_tmp.push_front(s_cur);
						_sdot_tmp.push_front(sdot_cur);
					}

					_s_tmp.pop_front();
					_sdot_tmp.pop_front();

					_s.insert(_s.end(), _s_tmp.begin(), _s_tmp.end());
					_sdot.insert(_sdot.end(), _sdot_tmp.begin(), _sdot_tmp.end());

					_s_tmp.clear();
					_sdot_tmp.clear();

					s_cur = _s.back();
					sdot_cur = _sdot.back();
					if (_switchPoint[_switchPoint.size() - 1]._id == SwitchPoint::SINGULAR)
					{
						// proceed following beta profile
						Real lambda = _switchPoint[_switchPoint.size() - 1]._lambda;
						for (int i = 0; i < numOfSPInt; i++)
						{
							s_cur += _ds;
							sdot_cur += lambda * _ds;
							_s.push_back(s_cur);
							_sdot.push_back(sdot_cur);
						}
					}
					beta_cur = determineAlphaBeta(s_cur, sdot_cur)(1);

					FI_SW = true;
					BI_SW = false;
				}
			}
		}

		// Step 3 : there exist two cases.
		s_cur = _sf-0.0001;
		sdot_cur = _vf / _dqds(_sf-0.0001).norm();

		_s_tmp.push_front(s_cur);
		_sdot_tmp.push_front(sdot_cur);

		if (!swiPoint_swi) // case 1. when can't find switching point until final s
		{
			LOGIF(_sdot.back() > sdot_cur, "step 3 error(!swiPoint_swi) : s_end has to be smaller than _sdot.back().");
			
			while (s_cur >= _s.back())
			{
				alphabeta = determineAlphaBeta(s_cur, sdot_cur);
				alpha_cur = alphabeta(0);
				backwardIntegrate(s_cur, sdot_cur, alpha_cur);
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
			_s.pop_back();
			_sdot.pop_back();
			alphabeta = determineAlphaBeta(s_cur, sdot_cur);
			alpha_cur = alphabeta(0);
			backwardIntegrate(s_cur, sdot_cur, alpha_cur);
			_s_tmp.push_front(s_cur);
			_sdot_tmp.push_front(sdot_cur);
		}

		_s_tmp.pop_front();
		_sdot_tmp.pop_front();

		_s.insert(_s.end(), _s_tmp.begin(), _s_tmp.end());
		_sdot.insert(_sdot.end(), _sdot_tmp.begin(), _sdot_tmp.end());
		_s_tmp.clear();
		_sdot_tmp.clear();

		calculateFinalTime();
		//calculateTorqueTrajectory();

		cout << "trajectory generation finished." << endl;
	}




	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////// HAVE TO DO //////////////////////////////////////////

	void TOPP::calculateTorqueTrajectory()
	{
		// TODO 수정중

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
		{
			q_data.col(i) = _q((*it_s));
			it_s++;
		}

		//std::vector<vector<Real>> torque_vec(_dof);
		//
		//for (int i = 0; i < _dof; i++)
		//{
		//	string torque_st = "C:/Users/crazy/Desktop/Time optimization/q";
		//	for (int j = 0; j < q_data.row(0).size(); j++)
		//		torque_vec[i].push_back(q_data(i, j));

		//	torque_st += to_string(i + 1);
		//	torque_st += ".txt";
		//	saveRealVector2txt(torque_vec[i], torque_st);
		//}
		//std::vector<Real> t_vec;
		//for (int i = 0; i < t.size(); i++)
		//	t_vec.push_back(t(i));
		//saveRealVector2txt(t_vec, "C:/Users/crazy/Desktop/Time optimization/t.txt");

		//cout << t << endl;
		//cout << "q_data.row(0) : " << q_data.row(0) << endl;
		//cout << "q_data.row(1) : " << q_data.row(1) << endl;
		//cout << "q_data column size : " << q_data.row(0).size() << endl;
		//cout << "t size : " << t.size() << endl;
		// make b-spline

		unsigned int degreeOfBSpline = 3;
		unsigned int orderOfBSpline = degreeOfBSpline + 1;
		unsigned int MaxNumOfCP = 15;
		unsigned int numData = _q_data.cols();
		//BSpline<-1, -1, -1> q = BSplineInterpolation(q_data, orderOfBSpline, 0, t(t.size()-1));
		//BSpline<-1, -1, -1> q = BSplineInterpolation(q_data, orderOfBSpline, t);
		BSpline<-1, -1, -1> q = BSplineFitting(q_data, orderOfBSpline, MaxNumOfCP, t(0), t(t.size() - 1));
		BSpline<-1, -1, -1> qdot = q.derivative();
		BSpline<-1, -1, -1> qddot = qdot.derivative();

		//cout << "q(0.05)" << endl;
		//cout << q(0.05) << endl;
		//cout << "q(0.10)" << endl;
		//cout << q(0.10) << endl;
		//cout << "q(0.15)" << endl;
		//cout << q(0.15) << endl;
		//cout << "q(0.20)" << endl;
		//cout << q(0.20) << endl;

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
}