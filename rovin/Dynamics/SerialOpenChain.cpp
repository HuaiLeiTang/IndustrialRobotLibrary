#include "SerialOpenChain.h"

namespace rovin {

	SerialOpenChain::SerialOpenChain() : _complete(false){}

	SerialOpenChain::~SerialOpenChain(){}

	LinkPtr SerialOpenChain::getLinkPtr(const unsigned int linkIdx)
	{
		LOGIF((linkIdx < _linkPtr.size()), "SerialOpenChain::getLinkPtr error : Link index is larger than number of Links");
		return _linkPtr[linkIdx];
	}

	MotorJointPtr SerialOpenChain::getMotorJointPtr(const unsigned int motorJointIdx)
	{
		LOGIF((motorJointIdx < _motorJointPtr.size()), "SerialOpenChain::getJointPtr error : Motorjoint index is larger than number of motorjoints");
		return _motorJointPtr[motorJointIdx];
	}

	bool SerialOpenChain::isComplete() const
	{
		return _complete;
	}

	const unsigned int SerialOpenChain::getNumOfLink() const
	{
		return _linkPtr.size();
	}

	const unsigned int SerialOpenChain::getNumOfJoint() const
	{
		return _motorJointPtr.size();
	}

	void SerialOpenChain::addLink(const LinkPtr & linkPtr)
	{
		_linkPtr.push_back(linkPtr);
		
	}

	void SerialOpenChain::deleteLink(unsigned int linkIdx)
	{
		LOGIF((linkIdx < _linkPtr.size()),"SerialOpenChain::deleteLink error : Link index is larger than number of Links");
		_linkPtr.erase(_linkPtr.begin() + linkIdx);
	}

	void SerialOpenChain::insertLink(unsigned int linkIdx, const LinkPtr& linkPtr)
	{
		LOGIF((linkIdx <= _linkPtr.size()), "SerialOpenChain::insertLink error : Link index is larger than number of Links");
		_linkPtr.insert(_linkPtr.begin() + linkIdx, linkPtr);
	}

	void SerialOpenChain::addMotorJoint(const MotorJointPtr & motorJointPtr)
	{
		_motorJointPtr.push_back(motorJointPtr);
	}

	void SerialOpenChain::deleteMotorJoint(unsigned int motorJointIdx)
	{
		LOGIF((motorJointIdx < _motorJointPtr.size()), "SerialOpenChain::deleteMotorJoint error : Motorjoint index is larger than number of motorjoints");
		_motorJointPtr.erase(_motorJointPtr.begin() + motorJointIdx);
	}

	void SerialOpenChain::insertMotorJoint(unsigned int motorJointIdx, const MotorJointPtr& motorJointPtr)
	{
		LOGIF((motorJointIdx <= _motorJointPtr.size()), "SerialOpenChain::insertMotorJoint error : Motorjoint index is larger than number of motorjoints");
		_motorJointPtr.insert(_motorJointPtr.begin() + motorJointIdx, motorJointPtr);
	}

	void SerialOpenChain::addMate(const unsigned int motorJointIdx, const SE3 & parentT, const SE3 & childT)
	{
		LOGIF((motorJointIdx < _motorJointPtr.size()), "SerialOpenChain::addMate error : Motorjoint index is larger than number of motorjoints");
		Mate mateInput(_motorJointPtr[motorJointIdx], parentT, childT);
		_Mate.push_back(mateInput);
	}

	void SerialOpenChain::completeAssembling()
	{
		unsigned int numOfmotorjoint = _motorJointPtr.size();
		unsigned int numOflink = _linkPtr.size();

		_socLink.resize(numOflink, SerialOpenChainLink());
		_socJoint.resize(numOfmotorjoint, SerialOpenChainJoint());

		LOGIF((numOfmotorjoint+1 == numOflink), "SerialOpenChain::completeAssembling error : numOfmotorjoint+1 != numOflink");

		// calculate each link frame, inertia and joint screw w.r.t base(fixed) frame
		SE3 T;
		se3 screw;
			
		T = _Mate[0].getParentT().inverse();
		screw = SE3::Ad(T) * (_Mate[0].getMotorJoint()->getAxis());
		_socJoint[0].setScrew(screw);
		
		T = T * _Mate[0].getChildT();
		_socLink[1].setM(T);

		Inertia inertia;
		inertia = _linkPtr[1]->getInertia();
		inertia.changeFrame(T);
		_socLink[1].setG(inertia);

		for (int i = 1; i < numOfmotorjoint; i++)
		{
			T = _socLink[i].getM() * _Mate[i].getParentT().inverse();
			screw = SE3::Ad(T) * (_Mate[i].getMotorJoint()->getAxis());
			_socJoint[i].setScrew(screw);

			T = T * _Mate[i].getChildT();
			_socLink[i + 1].setM(T);

			inertia = _linkPtr[i + 1]->getInertia();
			inertia.changeFrame(T);
			_socLink[i + 1].setG(inertia);
		}
		
		_complete = true;
	}

	StatePtr SerialOpenChain::makeState() const
	{
		LOGIF((0 < _motorJointPtr.size()), "SerialOpenChain::makeState() error : DOF should be larger than zero");
		StatePtr state(new State(_motorJointPtr.size()));

		// set baselink generalized velocity
		Vector6 V = Vector6::Zero();
		V[5] = 9.8;

		state->setLinkStateVel(0, V);

		return state;
	}

	/*!
	* \brief Serial Open Chain kinematics functions
	*/
	void SerialOpenChain::solveForwardKinematics(State & state)
	{
		LOGIF(_complete, "SerialOpenChain::solveForwardKinematics error : Assembly is not complete");

		if (state.checkStateInfoUpToDate(State::LINKS_POS))	{ return; }

		SE3 T;

		int jointsize = _motorJointPtr.size();

		if (state.checkStateInfoUpToDate(State::JOINTS_T_FROM_BASE))
		{
			for (int i = 0; i < jointsize; i++)
			{
				T = state.getJointStateAT(i) * _socLink[i + 1].getM();
				state.setLinkStateSE3(i + 1, T);
			}
		}
		else
		{
			for (int i = 0; i < jointsize; i++)
			{
				updateJointStateExponetial(state, i);
				if (i == 0)
					state.setJointStateAT(i, state.getJointStateT(i));
				else
				{
					T = state.getJointStateAT(i - 1) * state.getJointStateT(i);
					state.setJointStateAT(i, T);
				}
				T = state.getJointStateAT(i) * _socLink[i + 1].getM();
				state.setLinkStateSE3(i + 1, T);
			}
			state.updateStateInfoUpToDate(State::JOINTS_T_FROM_BASE, true);
		}
		state.updateStateInfoUpToDate(State::LINKS_POS, true);
	}

	void SerialOpenChain::solveDiffForwardKinematics(State & state)
	{
		LOGIF(_complete, "SerialOpenChain::solveForwardKinematics error : Assembly is not complete");

		if (state.checkStateInfoUpToDate(State::LINKS_VEL)) {
			return;
		}

		int jointsize = _motorJointPtr.size();
		se3 V;
		V = state.getLinkStateVel(0);
		
		if (state.checkStateInfoUpToDate(State::JOINTS_JACOBIAN))
		{
			// 이거 더 빠른지 check 하기..두 개의 계산 값이 같은지 확인하기
			for (int i = 0; i < jointsize; i++)
			{
				V += (state.getJointStateScrew(i) * state.getJointStateVel(i));
				state.setLinkStateVel(i + 1, V);
			}
		}
		else
		{
			for (int i = 0; i < jointsize; i++)
			{
				updateJointStateExponetial(state, i);
				V = SE3::Ad(state.getJointStateT(i).inverse(), V) + _socJoint[i].getScrew() * state.getJointStateVel(i);
				state.setLinkStateVel(i + 1, V);
			}
		}
		state.updateStateInfoUpToDate(State::LINKS_VEL, true);
	}

	void SerialOpenChain::solve2ndDiffForwardKinematics(State & state)
	{
		LOGIF(_complete, "SerialOpenChain::solve2ndDiffForwardKinematics error : Assembly is not complete");

		if (state.checkStateInfoUpToDate(State::LINKS_ACC)) {
			return;
		}

		int jointsize = _motorJointPtr.size();
		se3 V, VDot;
		V = state.getLinkStateVel(0);
		VDot = state.getLinkStateAcc(0);

		if (state.checkStateInfoUpToDate(State::JOINTS_JACOBIAN | State::JOINTS_JACOBIAN_DOT))
		{
			for (int i = 0; i < jointsize; i++)
			{
				VDot += (state.getJointStateScrewDot(i) * state.getJointStateVel(i) + state.getJointStateScrew(i) * state.getJointStateAcc(i));
				state.setLinkStateAcc(i + 1, VDot);
			}
		}
		else
		{
			solveDiffForwardKinematics(state);
			
			for (int i = 0; i < jointsize; i++)
			{
				updateJointStateExponetial(state, i);
				VDot = SE3::Ad(state.getJointStateT(i).inverse(), VDot) + SE3::ad(state.getLinkStateVel(i), _socJoint[i].getScrew()) * state.getJointStateVel(i)
					+ _socJoint[i].getScrew() * state.getJointStateAcc(i);
				state.setLinkStateAcc(i + 1, VDot);
			}
		}
		state.updateStateInfoUpToDate(State::LINKS_ACC, true);
	}

	void SerialOpenChain::solveJacobian(State & state)
	{
		LOGIF(_complete, "SerialOpenChain::solveJacobian error : Assembly is not complete");

		if (state.checkStateInfoUpToDate(State::JOINTS_JACOBIAN)) { return; }

		int jointsize = _motorJointPtr.size();
		SE3 T;
		for (int i = jointsize; i > 0; i--)
		{
			se3 scr = SE3::Ad(T, _socJoint[i - 1].getScrew());
			state.setJointStateScrew(i - 1, scr);
			updateJointStateExponetial(state, i - 1);
			T = state.getJointStateT(i - 1).inverse() * T;
		}
		state.updateStateInfoUpToDate(State::JOINTS_JACOBIAN, true);

	}

	void SerialOpenChain::solveJacobianDot(State & state)
	{
		LOGIF(_complete, "SerialOpenChain::solveJacobianDot error : Assembly is not complete");

		if (state.checkStateInfoUpToDate(State::JOINTS_JACOBIAN_DOT)) { return; }

		int jointsize = _motorJointPtr.size();
		se3 Jdot_i; /// (i-1)th column of Jdot
		se3 dJdq_qd; /// dJ(i-1)_dq(j-1)

		for (int i = jointsize; i > 0; i--) /// each column of Jdot
		{
			Jdot_i.setZero();
			for (int j = jointsize; j > i; j--) /// - dJ(i-1)_dq(j-1) * qdot
			{
				dJdq_qd = _socJoint[i - 1].getScrew();
				int k;
				for (k = jointsize; k > i; k--)
				{
					updateJointStateExponetial(state, k - 1);
					dJdq_qd = SE3::Ad(state.getJointStateT(k - 1).inverse(), dJdq_qd);
					if (k == (j - 1))
						dJdq_qd = SE3::ad(_socJoint[k].getScrew()) * dJdq_qd;;
				}
				Jdot_i -= dJdq_qd * state.getJointStateVel(j - 1);
			}
			state.setJointStateScrew(i - 1, Jdot_i);
		}
		state.updateStateInfoUpToDate(State::JOINTS_JACOBIAN_DOT, true);
	}

	void SerialOpenChain::updateJointStateExponetial(State & state, const unsigned int jointIndex)
	{
		if (!state.getJointState(jointIndex).checkJointInfoUpToDate(JointState::EXPOENTIAL))
		{
			state.setJointStateT(jointIndex, SE3::Exp(_socJoint[jointIndex].getScrew(), state.getJointStatePos(jointIndex)));
			state.getJointState(jointIndex).updateJointInfoUpToDate(JointState::EXPOENTIAL, true);
		}
	}


	// Dynamics

	void SerialOpenChain::solveInverDynamics(State & state, const dse3& endeffectorF)
	{


		// Forward Iteration
		solveForwardKinematics(state);
		solveDiffForwardKinematics(state);

		// Backward Iteration



	}

	void SerialOpenChain::solveFowardDynamics(State & state)
	{





	}



	// Mate class

	Mate::Mate(const MotorJointPtr & motorJoint, const SE3 & parentT, const SE3 & childT) : _parentT(parentT), _childT(childT)
	{
		_motorJoint = motorJoint->copy();
	}

	Mate::~Mate() {}

	void Mate::setMotorJoint(const MotorJointPtr & motorJoint)
	{
		_motorJoint = motorJoint->copy();
	}

	void Mate::setParentT(const SE3 & parentT)
	{
		_parentT = parentT;
	}

	void Mate::setChildT(const SE3 & childT)
	{
		_childT = childT;
	}

	const MotorJointPtr & Mate::getMotorJoint() const
	{
		return _motorJoint;
	}

	const SE3 & Mate::getParentT() const
	{
		return _parentT;
	}

	const SE3 & Mate::getChildT() const
	{
		return _childT;
	}

	// SerialOpenChainLink class

	SerialOpenChainLink::SerialOpenChainLink(const SE3 & M, const Inertia & G) : _M(M), _G(G) {}

	SerialOpenChainLink::~SerialOpenChainLink() {}

	void SerialOpenChainLink::setM(const SE3 & M)
	{
		_M = M;
	}

	void SerialOpenChainLink::setG(const Inertia & G)
	{
		_G = G;
	}

	const SE3 & SerialOpenChainLink::getM() const
	{
		return _M;
	}

	const Inertia & SerialOpenChainLink::getG() const
	{
		return _G;
	}

	SerialOpenChainJoint::SerialOpenChainJoint()
	{
		_screw.setZero();
	}

	// SerialOpenChainJoint class

	SerialOpenChainJoint::SerialOpenChainJoint(const se3 & screw) : _screw(screw)	{}

	SerialOpenChainJoint::~SerialOpenChainJoint()	{}

	void SerialOpenChainJoint::setScrew(const se3 & screw)
	{
		_screw = screw;
	}

	const se3 & SerialOpenChainJoint::getScrew() const
	{
		return _screw;
	}

}