#include "TOPP.h"

namespace rovin {

	TOPP::TOPP(const MatrixX & q_data, const SerialOpenChainPtr & soc, 
		const Real ds, const Real vi, const Real vf)
	{
		// _q_data dimesion : dof * number of data
		_q_data = q_data;
		_soc = soc;
		_ds = ds;
		_dof = _soc->getNumOfJoint();
		_vi = vi;
		_vf = vf;


		// insert torque constraint
		_torqueConstraint.resize(_dof * 2);
		for (int i = 0; i < _dof; i++)
		{
			_torqueConstraint[i] = -(_soc->getMotorJointPtr(i)->getLimitTorqueUpper());
			_torqueConstraint[i + _dof] = _soc->getMotorJointPtr(i)->getLimitAccLower();
		}

		// make spline
		_q = Spline(_q_data, 0, 1);
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

	VectorX& TOPP::calculateA(Real s, Real sdot)
	{
		VectorX a;
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

		a = tau - state->getJointStateTorque();
		return a;
	}

	bool TOPP::checkDynamicSingularity(const VectorX & a)
	{
		// TODO
		return false;
	}

	Vector2& TOPP::determineAlphaBeta(Real s, Real sdot)
	{
		VectorX a = calculateA(s, sdot);

		if (!checkDynamicSingularity(a))
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
			VectorX a_agg(_dof * 2);
			for (int i = 0; i < _dof; i++)
			{
				left_vec[i] = -tau[i] + _torqueConstraint[i];
				left_vec[i + _dof] = tau[i] - _torqueConstraint[i + _dof];
				a_agg[i] = a[i];
				a_agg[i + _dof] = a[i];
			}
			
			Vector2 result; // result[0] is alpha, result[1] is beta
			result[0] = -std::numeric_limits<Real>::max();
			result[1] = std::numeric_limits<Real>::max();

			for (int i = 1; i < _dof * 2; i++)
			{
				Real tmp = left_vec[i] / a_agg[i];
				if (a_agg[i] > 0) // upper bound beta
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

	Vector2 & TOPP::findNearestSwitchPoint(Real s, Real sdot)
	{
		// TODO: 여기에 반환 구문을 삽입합니다
		return Vector2();
	}

	void TOPP::generateTrajectory()
	{
		// s - sdot 을 어떤식으로 저장할 것인가..

		// initialization
		Real s_cur = 0;
		Real sdot_cur = _vi / _dqds(0).norm();
		_s.push_back(s_cur);
		_sdot.push_back(sdot_cur);

		Vector2& alphabeta = determineAlphaBeta(s_cur, sdot_cur);
		Real alpha_cur = alphabeta(0);
		Real beta_cur = alphabeta(1);
		
		// Step 1
		while (checkMVCCondition(alpha_cur, beta_cur))
		{

			
		}
		




	}

	



}