#include "TOPP.h"

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

		// insert torque constraint
		_torqueConstraint.resize(_dof * 2);
		for (int i = 0; i < _dof; i++)
		{
			_torqueConstraint[i] = -(_soc->getMotorJointPtr(i)->getLimitTorqueUpper());
			_torqueConstraint[i + _dof] = _soc->getMotorJointPtr(i)->getLimitAccLower();
		}

		// make b-spline
		_q = BSpline<-1, -1, -1>();
		_dqds = _q.derivative();
		_ddqdds = _dqds.derivative();

	}

	bool TOPP::checkMVCCondition(Real alpha, Real beta)
	{
		if (beta > alpha)
			return true;
		else
			return false;
	}

	bool TOPP::checkDynamicSingularity(const VectorX & a)
	{
		// TODO
		return false;
	}

	VectorX& TOPP::calculateA(Real s)
	{
		VectorX a(_dof * 2);
		VectorX q = _q(s);
		VectorX qs = _dqds(s);
		VectorX qdot(q.size());
		qdot.setZero();

		StatePtr state = _soc->makeState();
		state->setJointStatePos(q);
		state->setJointStateVel(qdot);
		state->setJointStateAcc(qs);

		_soc->solveInverseDynamics(*state);

		VectorX tau = state->getJointStateTorque();

		qs.setZero();

		state->setJointStateAcc(qs);

		_soc->solveInverseDynamics(*state);

		VectorX tmp_a = tau - state->getJointStateTorque();

		for (int i = 0; i < _dof; i++)
		{
			a[i] = tmp_a[i];
			a[i + _dof] = -tmp_a[i];
		}

		return a;
	}

	std::vector<VectorX>& TOPP::calculateBandC(Real s)
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
			b[i] = tmp_b[i];
			b[i + _dof] = -tmp_b[i];
			c[i] = tmp_c[i];
			c[i + _dof] = -tmp_c[i];
		}
		c = c + _torqueConstraint;

		result.push_back(b);
		result.push_back(c);

		return result;
	}

	Real TOPP::calculateMVCPoint(Real s)
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

	void TOPP::calculateFinalTime()
	{
		int integrateType = 1; // 1 : GQ, 2 : Euler

		if (integrateType == 1)
		{
			// Spline 만들고 다시 짜야징!
			// TODO
			GaussianQuadrature GQ(30, _si, _sf);

		}
		else if (integrateType == 2)
		{
			Real sum = 0;
			Real s_cur;
			while (_s.empty())
			{
				s_cur = _s.front();
				_s.pop_front();
				if (!_s.empty())
					sum += 1 / (_sdot.front()) * (_s.front() - s_cur);
				else
					sum += 1 / (_sdot.front()) * (_sf - s_cur);
				_sdot.front();
			}
			_tf_result = sum;
		}
		
	}

	void TOPP::calculateTorqueTrajectory()
	{
		// TODO
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

	Vector2& TOPP::determineAlphaBeta(Real s, Real sdot)
	{
		VectorX a = calculateA(s);

		// zero-inertia point check
		bool zero_inertia_swi = false;
		for (int i = 0; i < a.size(); i++)
		{
			if (a[i] < RealEps)
			{
				zero_inertia_swi = true;
				break;
			}
		}

		//if (!checkDynamicSingularity(a))
		if(!zero_inertia_swi)
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
			VectorX left_vec(_dof * 2);
			//VectorX a_agg(_dof * 2);
			for (int i = 0; i < _dof; i++)
			{
				left_vec[i] = -tau[i] + _torqueConstraint[i];
				left_vec[i + _dof] = tau[i] - _torqueConstraint[i + _dof];
				//a_agg[i] = a[i];
				//a_agg[i + _dof] = a[i];
			}
			
			Vector2 result; // result[0] is alpha, result[1] is beta
			result[0] = -std::numeric_limits<Real>::max();
			result[1] = std::numeric_limits<Real>::max();

			for (int i = 1; i < _dof * 2; i++)
			{
				Real tmp = left_vec[i] / a[i];
				if (a[i] > 0) // upper bound beta
				{
					if (tmp < result[1])
						result[1] = tmp;
				}
				else // lower bound alpha
				{
					if (tmp > result[0])
						result[0] = tmp;
				}
			}
			return result;
		}
		else // Dynamic  singularity
		{
			// TODO
			return Vector2();
		}
	}

	void TOPP::farwardIntegrate(Real & s, Real & sdot, Real sddot)
	{
		Real tmp = 2 * sddot*_ds + sdot*sdot;
		LOGIF((tmp > 0), "TOPP::farwardIntegrate error : the value has to be positive.");
		
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

		Real ds = 0.01;

		Real s_bef = s;
		Real sdot_bef = calculateMVCPoint(s_bef);
		Real s_cur = s + ds;
		Real sdot_cur = calculateMVCPoint(s_cur);
		Real s_next = s + 2 * ds;
		Real sdot_next = calculateMVCPoint(s_next);

		// 원래 topp 코드에는 tan_bef/tan_cur 로 discontinous 포인트 추가하는데 그거 이해 안감
		// 안넣어도 될것같은데 그럴거면 tan_bef/tan_cur 따로 구할 필요 없음
		//Real tan_bef = (sdot_cur - sdot_bef) / ds;
		//Real tan_cur = (sdot_next - sdot_cur) / ds;
		Real diff_bef = determineAlphaBeta(s_bef, sdot_bef)[0] / sdot_bef - (sdot_cur - sdot_bef) / ds;
		Real diff_cur = determineAlphaBeta(s_cur, sdot_cur)[0] / sdot_cur - (sdot_next - sdot_cur) / ds;

		VectorX a_bef = calculateA(s_bef);
		VectorX a_cur = calculateA(s_cur);
		std::vector<VectorX> bc_cur = calculateBandC(s_cur);


		// end criterion
		// TODO: check!
		if (s_next >= _sf)
			return false;


		while (true)
		{
			// 원래 코드는 singular -> tangent -> discontinuity 순서임... 왜지... 왜그런거지...
			// step 0: check whether a_k(s) is zero or not


			// step 1: discontinuous point check
			if (std::abs(sdot_next - sdot_cur) > 100 * std::abs(sdot_cur - sdot_bef))
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
			// step 1 -END-


			// step 2: singular point check
			// include calculating \lambda
			for (unsigned int i = 0; i < 2 * _dof; i++) // number of inequality
			{
				if (a_bef[i] * a_cur[i] <= 0)
				{
					// 딱 0이 아니면 원래 코드처럼 선형보간 해서 더 좋은 s를 찾는 과정이 필요할듯..
					// 여기 말고 위아래 step 1, 3에도 비슷한 알고리즘 필요할듯
					Real f = bc_cur[1][i] / bc_cur[0][i]; // \frac{c}{b}
					if (f < 0)
					{
						Real sdot_star = sqrt(-f);
						Real sdot_plus = calculateMVCPointExclude(s_cur, i);
						if (sdot_plus <= 0 && sdot_star < sdot_plus) // 첫번째 조건....
						{
							// calculate lambda
							Real diffeps = 1E-5;
							Real ap, bp, cp; // p for prime, differentiated to s
							ap = (calculateA(s_cur + diffeps)[i] - a_cur[i]) / diffeps;
							std::vector<VectorX> tmpbc = calculateBandC(s_cur + diffeps);
							bp = (tmpbc[0][i] - bc_cur[0][i]) / diffeps;
							cp = (tmpbc[1][i] - bc_cur[1][i]) / diffeps;

							Real lambda;

							if ((2 * bc_cur[0][i] + ap)*sdot_cur > 1E-10)
								lambda = (-bp*sdot_cur*sdot_cur - cp) / ((2 * bc_cur[0][i] + ap)*sdot_cur);
							else
								lambda = 0.0; // 이렇게 되면 어떡하징.. 원래 코드가 왜이렇지 흠 이거 들어오기는 할까..

							SwitchPoint sw(s_cur, sdot_cur, SwitchPoint::SINGULAR, lambda);
							_switchPoint.push_back(sw);
							return true;
						}
					}

				}
			}




			// step 3: tanget point check
			if ((diff_bef*diff_cur < 0) && (std::abs(diff_cur) < 1))
			{
				SwitchPoint sw(s_cur, sdot_cur, SwitchPoint::TANGENT, 0.0);
				_switchPoint.push_back(sw);
				return true;
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
		Real sdot_cur = _vi / _dqds(0).norm();
		_s.push_back(s_cur);
		_sdot.push_back(sdot_cur);

		Vector2& alphabeta = determineAlphaBeta(s_cur, sdot_cur);
		Real alpha_cur = alphabeta(0);
		Real beta_cur = alphabeta(1);

		// list used when backward integration
		std::list<Real> _s_tmp;
		std::list<Real> _sdot_tmp;

		bool FI_SW = true; ///< forward integration switch
		bool BI_SW = false; ///< backward integration switch
		bool I_SW = true; ///< integration switch

		bool swiPoint_swi;
		unsigned int numOfSPInt = 1;
		
		while (I_SW)
		{
			// Forward integration
			while (FI_SW)
			{
				// forward intergration
				farwardIntegrate(s_cur, sdot_cur, beta_cur); ///< update s_cur, sdot_cur

				// save trajectory points
				_s.push_back(s_cur);
				_sdot.push_back(sdot_cur);

				// calculate alpha and beta
				alphabeta = determineAlphaBeta(s_cur, sdot_cur);
				beta_cur = alphabeta(1);

				if (!checkMVCCondition(alpha_cur, beta_cur)) // case (a)
				{
					// fine nearest switch point
					swiPoint_swi = findNearestSwitchPoint(s_cur);

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
				else if (s_cur > _sf) ///< go to step 3, case (a), s_cur 가 무조건 s_end 넘어갔을 때!!
				{
					FI_SW = false;
					BI_SW = false;
					I_SW = false;
				}
				else if (sdot_cur < 1e-10) ///< Debugging 요소, case (a)
				{
					FI_SW = false;
					BI_SW = false;

				}
				
			}

			// Backward intergration
			while (BI_SW)
			{
				// calculate alpha and beta
				alphabeta = determineAlphaBeta(s_cur, sdot_cur);
				alpha_cur = alphabeta(0);

				// backward integration
				backwardIntegrate(s_cur, sdot_cur, alpha_cur); ///< update s_cur, sdot_cur

				// save trajectory points
				_s_tmp.push_front(s_cur);
				_sdot_tmp.push_front(sdot_cur);

				if (s_cur <= _s.back())
				{
					_s_tmp.pop_front();
					_sdot_tmp.pop_front();

					sdot_cur = (sdot_cur - _sdot_tmp.front()) / (s_cur - _s_tmp.front()) * (_s.back() - s_cur) + sdot_cur;
					s_cur = _s.back();

					LOGIF(_sdot.back() > sdot_cur,"Backward intergration error:_sdot.back() has to be larger than sdot_cur.");

					_s_tmp.push_front(s_cur);
					_sdot_tmp.push_front(sdot_cur);

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

					// switching point 저장 할 필요 없나요??

					_s.merge(_s_tmp);
					_sdot.merge(_sdot_tmp);

					_s_tmp.clear();
					_sdot_tmp.clear();

					FI_SW = true;
					BI_SW = false;
				}
			}
		}

		// Step 3 : there exist two cases.
		s_cur = _sf;
		sdot_cur = _vf / _dqds(_sf).norm();

		_s_tmp.push_front(s_cur);
		_sdot_tmp.push_front(sdot_cur);

		if (!swiPoint_swi) // case 1. when can't find switching point until final s
		{
			LOGIF(_sdot.back() > sdot_cur, "step 3 error(!swiPoint_swi) : s_end has to be smaller than _sdot.back().");
			
			while (s_cur <= _s.back())
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
		// switching point 저장 할 필요 없나요??
		_s.merge(_s_tmp);
		_sdot.merge(_sdot_tmp);
		_s_tmp.clear();
		_sdot_tmp.clear();
		_s.pop_back();
		_sdot.pop_back();

		// calculate tf and torque trajectory
		calculateFinalTime();
		calculateTorqueTrajectory();
	}
}