#include "GivenPathTimeOptimization.h"

using namespace std;

namespace rovin {
	GivenPathTimeOptimization::GivenPathTimeOptimization() {}

	GivenPathTimeOptimization::GivenPathTimeOptimization(const SerialOpenChainPtr& socRobotPtr, const std::vector<SE3, Eigen::aligned_allocator<SE3>>& givenPath)
	{
		_socRobotPtr = socRobotPtr;
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

	void GivenPathTimeOptimization::InputGivenPath(const std::vector<SE3, Eigen::aligned_allocator<SE3>>& givenPath)
	{
		_givenPath = givenPath;

		int pathsize = givenPath.size();
		_s.resize(pathsize);
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
		_q.resize(pathsize);
		for (int i = 0; i < pathsize; i++)
		{
			_socRobotPtr->solveInverseKinematics(*statePtr, _givenPath[i]);
			_q[i] = statePtr->getJointStatePos();
		}
		cout << "InversKin. Finish!" << endl;
		
		///< calculate qs: numerical differentiation
		_qs.resize(pathsize);
		for (int i = 0; i < pathsize; i++)
		{
			if (i == 0)
				_qs[i] = (_q[i + 1] - _q[i]) / (_s[i+1] - _s[i]);
			else if(i == (pathsize - 1))
				_qs[i] = (_q[i] - _q[i - 1]) / (_s[i] - _s[i - 1]);
			else
				_qs[i] = (_q[i + 1] - _q[i - 1]) / (_s[i + 1] - _s[i - 1]);
		}
		cout << "Prime73" << endl;

		///< calculate qss: numerical differentiation
		_qss.resize(pathsize);
		for (int i = 0; i < pathsize; i++)
		{
			if (i == 0)
				_qss[i] = (_qs[i + 1] - _qs[i]) / (_s[i + 1] - _s[i]);
			else if (i == (pathsize - 1))
				_qss[i] = (_qs[i] - _qs[i - 1]) / (_s[i] - _s[i - 1]);
			else
				_qss[i] = (_qs[i + 1] - _qs[i - 1]) / (_s[i + 1] - _s[i - 1]);
		}
		cout << "HYS love KWY" << endl;


	}

	void GivenPathTimeOptimization::MinMaxAcc(const unsigned int index, const Real& sdot, Real& min, Real& max)
	{

		StatePtr statePtr = _socRobotPtr->makeState();
		int dof = statePtr->getDof();

		statePtr->setJointStatePos(_q[index]);
		MatrixX M = _socRobotPtr->calculateMassMatrix(*statePtr);

		statePtr->setJointStateVel(_qs[index] * sdot);
		VectorX zerovector(dof);
		zerovector.setZero();
		statePtr->setJointStateAcc(zerovector);
		_socRobotPtr->solveInverseDynamics(*statePtr);

		VectorX h(dof);
		for (int i = 0; i < dof; i++)
			h(i) = statePtr->getJointStateTorque(i);

		VectorX c1(dof), c2(dof);
		c1 = M*_qs[index];
		c2 = M*_qss[index] * sdot*sdot + h;

		min = -RealMax;
		max = RealMax;
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