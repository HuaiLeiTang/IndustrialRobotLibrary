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

	VectorX GivenPathTimeOptimization::MinAcc(const Real & s, const Real & sdot)
	{

		return VectorX();
	}

	VectorX GivenPathTimeOptimization::MaxAcc(const Real & s, const Real & sdot)
	{
		return VectorX();
	}



}