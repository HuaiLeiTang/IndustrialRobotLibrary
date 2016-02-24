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

		// make spline
		_q = Spline(_q_data, _si, _sf);
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

	Real TOPP::calulateMVCPoint(Real s)
	{
		Real sdot_MVC = std::numeric_limits<Real>::max();

		Vector2 alphabeta = determineAlphaBeta(s, 0);
		if (alphabeta[0] > alphabeta[1])
			return 0;


		VectorX a = calculateA(s);
		std::vector<VectorX> BandC = calculateBandC(s);
		VectorX b = BandC[0];
		VectorX c = BandC[1];

		LOGIF(((a.size()==b.size()) && (a.size() == c.size()) && (b.size() == c.size())), "TOPP::calulateMVCPoint error : a, b, c vector size is wrong.")

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
		//
		// TODO: 여기에 반환 구문을 삽입합니다
		return true;
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

		bool FI_SW = true; ///< forward integration switch
		bool BI_SW = false; ///< backward integration switch
		bool I_SW = true; ///< integration switch

		bool swiPoint_swi;
		
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

				if (!checkMVCCondition(alpha_cur, beta_cur)) // go to backward integration
				{
					FI_SW = false;
					BI_SW = true;

					// fine nearest switch point
					swiPoint_swi = findNearestSwitchPoint(s_cur);

					s_cur = _switchPoint[_switchPoint.size() - 1]._s;
					sdot_cur = _switchPoint[_switchPoint.size() - 1]._sdot;

					// if swtich point doesn't exist until s_end --> go to step 3
					if (swiPoint_swi)
					{
						//std::cout << "Can not find switching point." << std::endl;
						FI_SW = false;
						BI_SW = false;
						I_SW = false;
					}
				}
				else if (sdot_cur < 1e-10) ///< Debugging 요소
				{
					FI_SW = false;
					BI_SW = false;

				}
				else if ((_sf - s_cur) < _ds) ///< go to step 3
				{
					FI_SW = false;
					BI_SW = false;
					I_SW = false;
				}
			}

			// Backward intergration
			while (BI_SW)
			{
				std::list<Real> _s_tmp;
				std::list<Real> _sdot_tmp;

				// save trajectory points
				_s_tmp.push_front(s_cur);
				_sdot_tmp.push_front(sdot_cur);

				// calculate alpha and beta
				alphabeta = determineAlphaBeta(s_cur, sdot_cur);
				alpha_cur = alphabeta(0);

				// backward integration
				backwardIntegrate(s_cur, sdot_cur, alpha_cur); ///< update s_cur, sdot_cur

				if (s_cur <= _s.back())
				{
					s_cur = _s.back();
				}
				

			}

		}

		// Step 3
		



	}
}