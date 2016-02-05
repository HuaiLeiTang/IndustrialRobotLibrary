#include "GivenPathTimeOptimization.h"

namespace rovin {
	GivenPathTimeOptimization::GivenPathTimeOptimization() {}

	GivenPathTimeOptimization::GivenPathTimeOptimization(const std::vector<SE3>& givenPath)
	{
		InputGivenPath(givenPath);
	}

	GivenPathTimeOptimization::~GivenPathTimeOptimization()	{}

	VectorX GivenPathTimeOptimization::getq(const unsigned int i) const
	{
		return _q[i];
	}

	VectorX GivenPathTimeOptimization::getqs(const unsigned int i) const
	{
		return _qs[i];
	}

	VectorX GivenPathTimeOptimization::getqss(const unsigned int i) const
	{
		return _qss[i];
	}

	SerialOpenChainPtr GivenPathTimeOptimization::getSOCPtr()
	{
		return _socRobotPtr;
	}

	void GivenPathTimeOptimization::InputGivenPath(const std::vector<SE3>& givenPath)
	{
		_givenPath = givenPath;

		int pathsize = givenPath.size();
		for (int i = 0; i < pathsize; i++)
			_s[i] = Real(i) / Real(pathsize - 1);

		solveInvKinAll();
	}

	void GivenPathTimeOptimization::solveMinimumTimeOptimization()
	{

	}

	void GivenPathTimeOptimization::solveInvKinAll()
	{
		int pathsize = _givenPath.size();
		LOGIF((pathsize != 0), "GivenPathTimeOptimization::solveInvKinAll error : No given path trajectory");
		
		StatePtr statePtr = _socRobotPtr->makeState();
		int dof = statePtr->getDof();

		///< calculate qs: inverse kinematics at each given location(SE3)
		for (int i = 0; i < pathsize; i++)
		{
			_socRobotPtr->solveInverseKinematics(*statePtr, _givenPath[i]);
			_q[i] = statePtr->getJointStatePos();
		}
		
		///< calculate qs: numerical differentiation
		for (int i = 0; i < pathsize; i++)
		{
			if (i == 0)
				_qs[i] = (_q[i + 1] - _q[i]) / (_s[i+1] - _s[i]);
			else if(i == (pathsize - 1))
				_qs[i] = (_q[i] - _q[i - 1]) / (_s[i] - _s[i - 1]);
			else
				_qs[i] = (_q[i + 1] - _q[i - 1]) / (_s[i + 1] - _s[i - 1]);
		}

		///< calculate qss numerical differentiation
		for (int i = 0; i < pathsize; i++)
		{
			if (i == 0)
				_qs[i] = (_q[i + 1] - _q[i]) / (_s[i + 1] - _s[i]);
			else if (i == (pathsize - 1))
				_qs[i] = (_q[i] - _q[i - 1]) / (_s[i] - _s[i - 1]);
			else
				_qs[i] = (_q[i + 1] - _q[i - 1]) / (_s[i + 1] - _s[i - 1]);
		}


	}

	void GivenPathTimeOptimization::MinMaxAcc(const unsigned int index, const Real& sdot, Real& min, Real& max)
	{
		int dof = _socRobotPtr->getNumOfJoint();

		State state(dof);
		state.setJointStatePos(_q[index]);
		MatrixX M(dof, dof);
		M = _socRobotPtr->calculateMassMatrix(state);

		state.setJointStateVel(_qs[index] * sdot);
		_socRobotPtr->solveInverseDynamics(state);
		VectorX h(dof);
		for (int i = 0; i < dof; i++)
			h(i) = state.getJointStateTorque(i);

		VectorX c1(dof), c2(dof);
		c1 = M*_qs[index];
		c2 = M*_qss[index] * sdot*sdot + h;

		min = 0;
		max = 0;
		Real mincur, maxcur;
		for (int i = 0; i < dof; i++)
		{
			if (c1[i] > 0)
			{
				mincur = (_socRobotPtr->getMotorJointPtr(i)->getLimitTorqueLower() - c2[i]) / c1[i];
				maxcur = (_socRobotPtr->getMotorJointPtr(i)->getLimitTorqueUpper() - c2[i]) / c1[i];
			}
			else
			{
				mincur = (_socRobotPtr->getMotorJointPtr(i)->getLimitTorqueUpper() - c2[i]) / c1[i];
				maxcur = (_socRobotPtr->getMotorJointPtr(i)->getLimitTorqueLower() - c2[i]) / c1[i];
			}

			if (mincur > min)
				min = mincur;
			if (maxcur < max)
				max = maxcur;
		}
	}


}