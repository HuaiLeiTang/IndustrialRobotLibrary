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
	
	void GivenPathTimeOptimization::setSerialOpenChainPtr(const SerialOpenChainPtr & socPtr)
	{
		_socRobotPtr = socPtr;
	}

	void GivenPathTimeOptimization::InputGivenPath(const std::vector<SE3, Eigen::aligned_allocator<SE3>>& givenPath)
	{
		_givenPath = givenPath;

		int pathsize = givenPath.size();
		_s.resize(pathsize);
		for (int i = 0; i < pathsize; i++)
			_s[i] = Real(i) / Real(pathsize - 1);

		_sdot.resize(pathsize);
		_sddot.resize(pathsize);

		solveInvKinAll();
	}

	void GivenPathTimeOptimization::solveMinimumTimeOptimization(const Real & sdot_i, const Real & sdot_f)
	{
		Real min, max;
		Real sdot;
		const int pathsize = _givenPath.size();  ///< pathsize

		const int InitStep = 0; ///< initial step for one forward iteration
		const int Finalstep = pathsize - 1;  ///< final step for one backward iteration

		int step; ///< 
		
		vector<Real> sdotarray_fwd;
		vector<Real> sdotarray_bwd;
		sdotarray_fwd.resize(pathsize);
		sdotarray_bwd.resize(pathsize);
		for (int i = 0; i < pathsize; i++)
		{
			sdotarray_fwd[i] = RealMax;
			sdotarray_bwd[i] = -RealMax;
		}

		int binary_num;

		int fw_startstep = InitStep;
		int bw_startstep;

		bool isfeasible;
		int SwitchingPoint1, SwitchingPoint2;
		Real max_minus_min;
		Real minimumMinMaxDiff;

		Real fw_startsdot = sdot_i;
		while (true)
		{
			///< forward search
			step = fw_startstep;
			sdot = fw_startsdot;
			while (true)
			{
				sdotarray_fwd[step] = sdot;
				MinMaxAcc(step, sdot, min, max);

				if (min > max)
				{
					binary_num = pathsize - step;
					break;
				}

				sdot = sqrt(sdot*sdot + 2 * max * (_s[step + 1] - _s[step]));
				step++;

				if (step == Finalstep)
				{
					sdotarray_fwd[step] = sdot;
					break;
				}
			}


			///< backward search iteration
			bw_startstep = Finalstep;
			while(binary_num != 0)
			{
				step = bw_startstep;
				sdot = sdot_f;
				if (binary_num > 1)
					binary_num = (binary_num + 1) / 2;
				else
					binary_num = 0;

				minimumMinMaxDiff = RealMax;
				while (true) ///< one backward search
				{
					sdotarray_bwd[step] = sdot;
					MinMaxAcc(step, sdot, min, max);
					max_minus_min = max - min;

					if ((max_minus_min < minimumMinMaxDiff) && (max_minus_min > 0))
					{
						minimumMinMaxDiff = max_minus_min;
						SwitchingPoint2 = step;
					}

					if (sdot > sdotarray_fwd[step]) ///< feasible.
					{
						bw_startstep += binary_num;
						isfeasible = true;
						SwitchingPoint1 = step;
						break;
					}

					if (max_minus_min < 0) ///< not feasible.
					{
						bw_startstep -= binary_num;
						isfeasible = false;
						break;
					}


					sdot = sqrt(sdot*sdot - 2 * min * (_s[step] - _s[step - 1]));
					step--;

				} ///< one backward search end
				
				if (bw_startstep > Finalstep) ///< if first backward search is feasible.
				{
					SwitchingPoint2 = Finalstep;
					break;
				}

			}///< backward search iteration end


			///< repeat when last backward search is not feasible.
			if (!isfeasible)
			{
				step = bw_startstep - 1;
				sdot = sdot_f;
				minimumMinMaxDiff = RealMax;
				while (true)
				{
					sdotarray_bwd[step] = sdot;
					MinMaxAcc(step, sdot, min, max);
					max_minus_min = max - min;

					LOGIF((max_minus_min > 0), "something is wrong");
					if (max_minus_min <= 0)
					{
						cout << bw_startstep << endl;
						cout << step << endl;
						cout << binary_num << endl;
					}

					if ((minimumMinMaxDiff > max_minus_min) && (max_minus_min > 0))
					{
						minimumMinMaxDiff = max_minus_min;
						SwitchingPoint2 = step;
					}

					if (sdot > sdotarray_fwd[step]) ///< feasible.
					{
						SwitchingPoint1 = step;
						break;
					}

					sdot = sqrt(sdot*sdot - 2 * min * (_s[step] - _s[step - 1]));
					step--;
				}
				cout << "backward once again" << endl;
			}


			cout << "fw_startstep: " << fw_startstep << "  SwitchingPoint1: " << SwitchingPoint1 << "  SwitchingPoint2: " << SwitchingPoint2 << endl;
			//cout << "  minimumMinMaxDiff: " << minimumMinMaxDiff << endl;

			///< Save one pair of forward-backward trajectory.
			for (int i = fw_startstep; i < SwitchingPoint1 + 1; i++)
				_sdot[i] = sdotarray_fwd[i];

			for (int i = SwitchingPoint1 + 1; i < SwitchingPoint2 + 1; i++)
				_sdot[i] = sdotarray_bwd[i];

			for (int i = fw_startstep; i < SwitchingPoint2 + 1; i++)
				cout << "i: " << i << "  s: " << _s[i] << "   sdot: " << _sdot[i] << endl;

			fw_startstep = SwitchingPoint2;
			fw_startsdot = _sdot[fw_startstep];

			if (binary_num != 0) ///< All sdot search is done.
				break;

		}

		cout << "End!!!!!!!!!!" << endl;
		for (int i = 0; i < pathsize; i++)
		{
			cout << "i: " << i << "  s: " << _s[i] << "   sdot: "<< _sdot[i] << endl;
		}
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

		VectorX h = statePtr->getJointStateTorque();

		VectorX c1(dof), c2(dof);
		c1 = M*_qs[index];
		c2 = M*_qss[index] * sdot*sdot + h;
		
		min = -RealMax;
		max = RealMax;
		Real mincur, maxcur;

		for (int i = 0; i < dof; i++)
		{
			if (c1(i) > 0)
			{
				mincur = (_socRobotPtr->getMotorJointPtr(i)->getLimitTorqueLower() - c2(i)) / c1(i);
				maxcur = (_socRobotPtr->getMotorJointPtr(i)->getLimitTorqueUpper() - c2(i)) / c1[i];
			}
			else
			{
				mincur = (_socRobotPtr->getMotorJointPtr(i)->getLimitTorqueUpper() - c2(i)) / c1(i);
				maxcur = (_socRobotPtr->getMotorJointPtr(i)->getLimitTorqueLower() - c2(i)) / c1(i);
			}

			if (mincur > min)
				min = mincur;
			if (maxcur < max)
				max = maxcur;
		}

	}


}