#include "GivenPathTimeOptimization.h"

namespace rovin {
	GivenPathTimeOptimization::GivenPathTimeOptimization() {}

	GivenPathTimeOptimization::GivenPathTimeOptimization(const std::vector<SE3>& givenPath)
	{
		setGivenPath(givenPath);
	}

	GivenPathTimeOptimization::~GivenPathTimeOptimization()	{}

	VectorX GivenPathTimeOptimization::getq(const unsigned int i)
	{
		return _q[i];
	}

	VectorX GivenPathTimeOptimization::getqs(const unsigned int i)
	{
		return _qs[i];
	}

	VectorX GivenPathTimeOptimization::getqss(const unsigned int i)
	{
		return _qss[i];
	}

	void GivenPathTimeOptimization::setGivenPath(const std::vector<SE3>& givenPath)
	{
		_givenPath = givenPath;

		int pointsize = givenPath.size();
		for (int i = 0; i < pointsize; i++)
			_s[i] = i / Real(pointsize - 1);
	}

	void GivenPathTimeOptimization::solveMinimumTimeOptimization()
	{

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