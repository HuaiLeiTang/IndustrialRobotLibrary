#include "SerialOpenChain.h"

using namespace std;

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

	const LinkPtr& SerialOpenChain::getLinkPtr(const unsigned int linkIdx) const
	{
		LOGIF((linkIdx < _linkPtr.size()), "SerialOpenChain::getLinkPtr error : Link index is larger than number of Links");
		return _linkPtr[linkIdx];
	}

	const MotorJointPtr& SerialOpenChain::getMotorJointPtr(const unsigned int motorJointIdx) const
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

		for (unsigned int i = 1; i < numOfmotorjoint; i++)
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
		Vector6 VDot = Vector6::Zero();
		VDot[5] = 9.8;

		state->setLinkStateAcc(0, VDot);

		return state;
	}

	/*!
	* \brief Serial Open Chain kinematics functions
	*/
	void SerialOpenChain::solveForwardKinematics(State & state)
	{
		LOGIF(_complete, "SerialOpenChain::solveForwardKinematics error : Assembly is not complete");

		if (state.checkStateInfoUpToDate(State::LINKS_POS))	{ return; }

		updateAccumulatedT(state);

		SE3 T;
		int jointsize = _motorJointPtr.size();

		for (int i = 0; i < jointsize; i++)
		{
			T = state.getJointStateAT(i) * _socLink[i + 1].getM();
			state.setLinkStateSE3(i + 1, T);
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
				VDot = SE3::Ad(state.getJointStateT(i).inverse(), VDot) + SE3::ad(state.getLinkStateVel(i + 1), _socJoint[i].getScrew()) * state.getJointStateVel(i)
					+ _socJoint[i].getScrew() * state.getJointStateAcc(i);
				state.setLinkStateAcc(i + 1, VDot);
			}
		}
		state.updateStateInfoUpToDate(State::LINKS_ACC, true);
	}

	void SerialOpenChain::solveInverseKinematics(State & state, const SE3& goalT)
	{
		int index_endeffector = _linkPtr.size() - 1;
		int jointsize = _motorJointPtr.size();

		se3 S;
		MatrixX J(6, jointsize);
		SE3 goalTmod_inv;
		VectorX delta;

		goalTmod_inv = _socLink[index_endeffector].getM() * goalT.inverse();

		while (true)
		{
			updateAccumulatedT(state);

			if((S = SE3::Log(goalTmod_inv * state.getJointStateAT(index_endeffector - 1))).norm() < InverseKinematicsExitCondition)
				break;

			solveJacobian(state);
			for (int i = 0; i < jointsize; i++)
				J.col(i) = state.getJointStateScrew(i);

			delta = pInv(J) * (-S);

			for (int i = 0; i < jointsize; i++)
				state.addJointStatePos(i, delta[i]);
		}
	}

	void SerialOpenChain::solveJacobian(State & state)
	{
		LOGIF(_complete, "SerialOpenChain::solveJacobian error : Assembly is not complete");

		if (state.checkStateInfoUpToDate(State::JOINTS_JACOBIAN)) { return; }

		int jointsize = _motorJointPtr.size();
		SE3 T;
		for (int i = jointsize; i > 0; i--)
		{
			se3 scr = SE3::Ad(T.inverse(), _socJoint[i - 1].getScrew());
			state.setJointStateScrew(i - 1, scr);
			updateJointStateExponetial(state, i - 1);
			T = state.getJointStateT(i - 1) * T;
		}
		state.updateStateInfoUpToDate(State::JOINTS_JACOBIAN, true);

	}

	void SerialOpenChain::solveJacobianDot(State & state)
	{
		LOGIF(_complete, "SerialOpenChain::solveJacobianDot error : Assembly is not complete");

		if (state.checkStateInfoUpToDate(State::JOINTS_JACOBIAN_DOT)) { return; }

		int jointsize = _motorJointPtr.size();
		Matrix6 ad_sum = Matrix6::Zero();

		for (int i = jointsize; i > 0; i--) /// each column of Jdot
		{
			ad_sum += SE3::ad(state.getJointStateScrew(i - 1)) * state.getJointStateVel(i - 1);
			state.setJointStateScrewDot(i - 1, -ad_sum * state.getJointStateScrew(i - 1));
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

	void SerialOpenChain::updateAccumulatedT(State & state)
	{

		if (!state.checkStateInfoUpToDate(State::JOINTS_T_FROM_BASE))
		{
			int jointsize = _motorJointPtr.size();
			SE3 T;

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
			}

			state.updateStateInfoUpToDate(State::JOINTS_T_FROM_BASE, true);
		}
	}


	// Dynamics

	void SerialOpenChain::solveInverseDynamics(State & state, const dse3& endeffectorF)
	{
		// Forward Iteration
		solveForwardKinematics(state);
		solveDiffForwardKinematics(state);
		solve2ndDiffForwardKinematics(state);

		// Backward Iteration
		int jointsize = _motorJointPtr.size();
		dse3 F_end = (SE3::Ad(_socLink[_motorJointPtr.size()].getM().inverse())).transpose() * endeffectorF; ///< transformation end-effector dse3
		dse3 F;

		for (int i = jointsize; i > 0; i--)
		{
			const Matrix6& G = static_cast<const Matrix6&>(_socLink[i].getG());
			F = G * state.getLinkStateAcc(i) - SE3::adTranspose(state.getLinkStateVel(i), G * state.getLinkStateVel(i));
			if (i == jointsize)
				F += F_end;
			else
				F += SE3::Ad(state.getJointStateT(i).inverse()).transpose() * state.getJointStateConstraintF(i);
			state.setJointStateConstraintF(i - 1, F);

			// torque
			Real tau = F.dot(_socJoint[i - 1].getScrew());
			tau += state.getJointStatePos(i - 1)*_motorJointPtr[i - 1]->getSpringConstant() + state.getJointStateVel(i - 1)*_motorJointPtr[i - 1]->getDamperConstant();
			if (RealBigger(state.getJointStateVel(i - 1), 0))
				tau += _motorJointPtr[i - 1]->getCoulombFrictionConstant();
			else if (RealLess(state.getJointStateVel(i - 1), 0))
				tau -= _motorJointPtr[i - 1]->getCoulombFrictionConstant();
			state.setJointStateTorque(i-1, tau);
		}
	}

	void SerialOpenChain::solveFowardDynamics(State & state)
	{





	}

	MatrixX SerialOpenChain::differentiateInverseDynamics(State & state, const MatrixX & dqdp, const MatrixX & dqdotdp, const MatrixX & dqddotdp)
	{
		int dof = dqdp.rows(); ///< Robot degree of freedom
		int pN = dqdp.cols(); ///< number of parameters
		int linkN = getNumOfLink();

		// variables for forward iteration
		std::vector<MatrixX> dVdp(linkN), dVdotdp(linkN);

		// variables for backward iteration
		std::vector<MatrixX> dFdp(dof);
		MatrixX dtaudp(dof, pN);

		// Initialization
		dVdp[0] = MatrixX::Zero(6, pN);
		dVdotdp[0] = MatrixX::Zero(6, pN);
		dFdp[0] = MatrixX::Zero(6, pN);

		// Forward recursion
		MatrixX curdVdp = MatrixX::Zero(6, pN);
		MatrixX currdVdotdp = MatrixX::Zero(6, pN);

		for (int i = 0; i < dof; i++)
		{
			// calculate current dVdp, dimension 6 * pN
			curdVdp = SE3::Ad(state.getJointStateT(i).inverse()) * curdVdp + 
				_socJoint[i].getScrew() * dqdotdp.row(i) -
				SE3::ad(_socJoint[i].getScrew(), state.getLinkStateVel(i + 1)) * dqdp.row(i);

			// calculate current dVdotdq, dimension 6 * pN
			currdVdotdp = SE3::Ad(state.getJointStateT(i).inverse()) * currdVdotdp +
				_socJoint[i].getScrew() * dqddotdp.row(i) -
				SE3::ad(_socJoint[i].getScrew(), state.getLinkStateVel(i + 1)) * dqdotdp.row(i) -
				SE3::ad(_socJoint[i].getScrew()) * curdVdp * state.getJointStateVel(i) -
				SE3::ad(_socJoint[i].getScrew()) * SE3::Ad(state.getJointStateT(i).inverse(), state.getLinkStateAcc(i)) * dqdp.row(i);

			// save data
			dVdp[i + 1] = curdVdp;
			dVdotdp[i + 1] = currdVdotdp;
		}

		// Backward recursion
		MatrixX curdFdp = MatrixX::Zero(6, pN);

		for (int i = dof; i > 0; i--)
		{
			// calculate current dFdp, dimension 6 * pN
			const Matrix6& G = static_cast<const Matrix6&>(_socLink[i].getG());
			curdFdp += G * dVdotdp[i] -
				SE3::adTranspose(state.getLinkStateVel(i)) * G * dVdp[i];
			se3 tmp = G * state.getLinkStateVel(i);
			for (int k = 0; k < pN; k++)
			{
				curdFdp.col(k) -= SE3::adTranspose(dVdp[i].col(k), tmp);
			}

			// calculate current dtaudp, dimension 1 * pN
			se3 screw = _socJoint[i - 1].getScrew();
			dtaudp.row(i - 1) = screw.transpose() * curdFdp + _motorJointPtr[i - 1]->getSpringConstant()*dqdp.row(i - 1) + _motorJointPtr[i - 1]->getDamperConstant()*dqdotdp.row(i - 1);

			curdFdp = SE3::Ad(state.getJointStateT(i - 1).inverse()).transpose() * (SE3::adTranspose(-(_socJoint[i - 1].getScrew()), state.getJointStateConstraintF(i - 1)) * dqdp.row(i - 1) + curdFdp);
		}

		return dtaudp;
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