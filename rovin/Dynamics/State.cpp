#include "State.h"

namespace rovin {

	//////////////////////////////////////////////////////////////
	//						 STATE CLASS			     		//
	//////////////////////////////////////////////////////////////

	State::State(unsigned int dof) : _dof(dof), _stateInfoUpToDate(0)
	{
		for (int i = 0; i < _dof+1; i++)
		{
			_linkState.push_back(LinkState());
		}

		for (int i = 0; i < _dof; i++)
		{
			_jointState.push_back(JointState());
		}
	}

	State::~State() {}

	const int State::getDof() const
	{
		return _dof;
	}

	JointState & State::getJointState(const unsigned int jointIndex)
	{
		return _jointState[jointIndex];
	}

	LinkState & State::getLinkState(const unsigned int linkIndex)
	{
		return _linkState[linkIndex];
	}

	const JointState & State::getJointState(const unsigned int jointIndex) const
	{
		return _jointState[jointIndex];
	}

	const LinkState & State::getLinkState(const unsigned int linkIndex) const
	{
		return _linkState[linkIndex];
	}

	const Real State::getJointStatePos(const unsigned int jointIndex) const
	{
		return _jointState[jointIndex].getJointPos();
	}

	const Real State::getJointStateVel(const unsigned int jointIndex) const
	{
		return _jointState[jointIndex].getJointVel();
	}

	const Real State::getJointStateAcc(const unsigned int jointIndex) const
	{
		return _jointState[jointIndex].getJointAcc();
	}

	const VectorX State::getJointStatePos() const
	{
		VectorX q(_dof);
		for (unsigned int i = 0; i < _dof; i++)
		{
			q(i) = _jointState[i].getJointPos();
		}
		return q;
	}

	const VectorX State::getJointStateVel() const
	{
		VectorX qdot(_dof);
		for (unsigned int i = 0; i < _dof; i++)
		{
			qdot(i) = _jointState[i].getJointVel();
		}
		return qdot;
	}

	const VectorX State::getJointStateAcc() const
	{
		VectorX qddot(_dof);
		for (unsigned int i = 0; i < _dof; i++)
		{
			qddot(i) = _jointState[i].getJointAcc();
		}
		return qddot;
	}

	const se3 & State::getJointStateScrew(const unsigned int jointIndex) const
	{
		return _jointState[jointIndex].getJointScrew();
	}

	const se3 & State::getJointStateScrewDot(const unsigned int jointIndex) const
	{
		return _jointState[jointIndex].getJointScrewDot();
	}

	const SE3& State::getJointStateT(const unsigned int jointIndex) const
	{
		return _jointState[jointIndex].getJointT();
	}

	const SE3& State::getJointStateAT(const unsigned int jointIndex) const
	{
		return _jointState[jointIndex].getJointAccumulatedT();
	}

	const dse3 & State::getJointStateConstraintF(const unsigned int jointIndex) const
	{
		return _jointState[jointIndex].getJointConstraintF();
	}

	const Real State::getJointStateTorque(const unsigned int jointIndex) const
	{
		return _jointState[jointIndex].getJointTorque();
	}

	const VectorX State::getJointStateTorque() const
	{
		VectorX tau(_dof);
		for (unsigned int i = 0; i < _dof; i++)
		{
			tau(i) = _jointState[i].getJointTorque();
		}
		return tau;
	}

	const SE3& State::getLinkStateSE3(const unsigned int linkIndex) const
	{
		return _linkState[linkIndex].getLinkSE3();
	}

	const se3 & State::getLinkStateVel(const unsigned int linkIndex) const
	{
		return _linkState[linkIndex].getLinkVel();
	}

	const se3 & State::getLinkStateAcc(const unsigned int linkIndex) const
	{
		return _linkState[linkIndex].getLinkVelDot();
	}

	void State::setJointStatePos(JointState & jointstate, const Real q)
	{
		jointstate.setJointPos(q);
		updateStateInfoUpToDate(LINKS_POS | LINKS_VEL | LINKS_ACC | JOINTS_T_FROM_BASE | JOINTS_JACOBIAN | JOINTS_JACOBIAN_DOT, false);
	}

	// joint 값 vector 로 받는 함수 추가

	void State::setJointStatePos(const unsigned int jointIndex, const Real q)
	{
		setJointStatePos(_jointState[jointIndex], q);
	}

	void State::setJointStateVel(JointState & jointstate, const Real qdot)
	{
		jointstate.setJointVel(qdot);
		updateStateInfoUpToDate(JOINTS_JACOBIAN_DOT | LINKS_VEL | LINKS_ACC, false);
	}

	void State::setJointStateVel(const unsigned int jointIndex, const Real qdot)
	{
		setJointStateVel(_jointState[jointIndex], qdot);
	}

	void State::setJointStateAcc(JointState & jointstate, const Real qddot)
	{
		jointstate.setJointAcc(qddot);
		updateStateInfoUpToDate(LINKS_ACC, false);
	}

	void State::setJointStateAcc(const unsigned int jointIndex, const Real qddot)
	{
		setJointStateAcc(_jointState[jointIndex], qddot);
	}

	void State::setJointStatePos(const VectorX& q)
	{
		for (unsigned int i = 0; i < _dof; i++)
		{
			setJointStatePos(_jointState[i], q(i));
		}
	}

	void State::setJointStateVel(const VectorX& qdot)
	{
		for (unsigned int i = 0; i < _dof; i++)
		{
			setJointStateVel(_jointState[i], qdot(i));
		}
	}

	void State::setJointStateAcc(const VectorX& qddot)
	{
		for (unsigned int i = 0; i < _dof; i++)
		{
			setJointStateAcc(_jointState[i], qddot(i));
		}
	}

	void State::setJointStateScrew(const unsigned int jointIndex, const se3 & screw)
	{
		_jointState[jointIndex].setJointScrew(screw);
	}

	void State::setJointStateScrewDot(const unsigned int jointIndex, const se3 & screwDot)
	{
		_jointState[jointIndex].setJointScrewDot(screwDot);
	}

	void State::setJointStateT(const unsigned int jointIndex, const SE3 & T)
	{
		_jointState[jointIndex].setJointT(T);
	}

	void State::setJointStateAT(const unsigned int jointIndex, const SE3 & T)
	{
		_jointState[jointIndex].setJointAccumulatedT(T);
	}

	void State::setJointStateConstraintF(const unsigned int jointIndex, const dse3 & constraintF)
	{
		_jointState[jointIndex].setJointConstraintF(constraintF);
	}

	void State::setJointStateTorque(const unsigned int jointIndex, const Real tau)
	{
		_jointState[jointIndex].setJointTorque(tau);
	}

	void State::setJointStateTorque(const VectorX& tau)
	{
		for (unsigned int i = 0; i < _dof; i++)
		{
			setJointStateTorque(i, tau(i));
		}
	}

	void State::setLinkStateSE3(const unsigned int linkIndex, const SE3 & T)
	{
		_linkState[linkIndex].setLinkSE3(T);
	}

	void State::setLinkStateVel(const unsigned int linkIndex, const se3 & V)
	{
		_linkState[linkIndex].setLinkVel(V);
	}

	void State::setLinkStateAcc(const unsigned int linkIndex, const se3 & VDot)
	{
		_linkState[linkIndex].setLinkVelDot(VDot);
	}

	void State::addJointStatePos(JointState & jointstate, const Real q)
	{
		jointstate.addJointPos(q);
		updateStateInfoUpToDate(LINKS_POS | LINKS_VEL | LINKS_ACC | JOINTS_T_FROM_BASE | JOINTS_JACOBIAN | JOINTS_JACOBIAN_DOT, false);
	}

	void State::addJointStatePos(const unsigned int jointIndex, const Real q)
	{
		addJointStatePos(_jointState[jointIndex], q);
	}

	void State::addJointStateVel(JointState & jointstate, const Real qdot)
	{
		jointstate.addJointVel(qdot);
		updateStateInfoUpToDate(JOINTS_JACOBIAN_DOT | LINKS_VEL | LINKS_ACC, false);
	}

	void State::addJointStateVel(const unsigned int jointIndex, const Real qdot)
	{
		addJointStateVel(_jointState[jointIndex], qdot);
	}

	void State::addJointStateAcc(JointState & jointstate, const Real qddot)
	{
		jointstate.addJointAcc(qddot);
		updateStateInfoUpToDate(LINKS_ACC, false);
	}

	void State::addJointStateAcc(const unsigned int jointIndex, const Real qddot)
	{
		addJointStateAcc(_jointState[jointIndex], qddot);
	}

	void State::addJointStatePos(const VectorX& q)
	{
		for (unsigned int i = 0; i < _dof; i++)
		{
			addJointStatePos(_jointState[i], q(i));
		}
	}

	void State::addJointStateVel(const VectorX& qdot)
	{
		for (unsigned int i = 0; i < _dof; i++)
		{
			addJointStateVel(_jointState[i], qdot(i));
		}
	}

	void State::addJointStateAcc(const VectorX& qddot)
	{
		for (unsigned int i = 0; i < _dof; i++)
		{
			addJointStateAcc(_jointState[i], qddot(i));
		}
	}

	bool State::checkStateInfoUpToDate(int infoIdx)
	{
		return infoIdx == (infoIdx & _stateInfoUpToDate);
	}

	void State::updateStateInfoUpToDate(int infoIdx, bool upToDate)
	{
		if (upToDate)
			_stateInfoUpToDate |= infoIdx;
		else
			_stateInfoUpToDate &= ~infoIdx;
	}
	
	//////////////////////////////////////////////////////////////
	//						LINK STATE CLASS					//
	//////////////////////////////////////////////////////////////

	LinkState::LinkState()
	{
		_V.setZero();
		_VDot.setZero();
		_Ja.setZero();
		_b.setZero();
	}

	LinkState::~LinkState() {}

	void LinkState::setLinkSE3(const SE3 & T)
	{
		_T = T;
	}
	void LinkState::setLinkVel(const se3 & V)
	{
		_V = V;
	}
	void LinkState::setLinkVelDot(const se3 & VDot)
	{
		_VDot = VDot;
	}
	void LinkState::setLinkArtInertia(const Matrix6 & Ja)
	{
		_Ja = Ja;
	}
	void LinkState::setLinkBiasforce(const dse3 & b)
	{
		_b = b;
	}
	const SE3 & LinkState::getLinkSE3() const
	{
		return _T;
	}
	const se3 & LinkState::getLinkVel() const
	{
		return _V;
	}
	const se3 & LinkState::getLinkVelDot() const
	{
		return _VDot;
	}
	const Matrix6 & LinkState::getLinkArtInertia() const
	{
		return _Ja;
	}
	const dse3 & LinkState::getLinkBiasforce() const
	{
		return _b;
	}

	LinkStatePtr LinkState::copy() const
	{
		LinkStatePtr clone(new LinkState());
		clone->setLinkSE3(_T);
		clone->setLinkVel(_V);
		clone->setLinkVelDot(_VDot);
		clone->setLinkArtInertia(_Ja);
		clone->setLinkBiasforce(_b);
		return clone;
	}


	//////////////////////////////////////////////////////////////
	//						JOINT STATE CLASS					//
	//////////////////////////////////////////////////////////////

	JointState::JointState() : _tau(0), _q(0), _qdot(0), _qddot(0), _jointInfoUpToDate(0)
	{
		_screw.setZero();
		_screwDot.setZero();
		_constraintF.setZero();
	}

	JointState::~JointState() {}

	void JointState::setJointTorque(const Real tau)
	{
		_tau = tau;
	}

	void JointState::setJointT(const SE3 & T)
	{
		_T = T;
	}

	void JointState::setJointAccumulatedT(const SE3 & accumulatedT)
	{
		_accumulatedT = accumulatedT;
	}

	void JointState::setJointScrew(const se3 & screw)
	{
		_screw = screw;
	}

	void JointState::setJointScrewDot(const se3 & screwdot)
	{
		_screwDot = screwdot;
	}

	void JointState::setJointConstraintF(const dse3 & constraintF)
	{
		_constraintF = constraintF;
	}

	void JointState::setJointPos(const Real q)
	{
		_q = q;
		updateJointInfoUpToDate(EXPOENTIAL, false);
	}

	void JointState::setJointVel(const Real qdot)
	{
		_qdot = qdot;
	}

	void JointState::setJointAcc(const Real qddot)
	{
		_qddot = qddot;
	}

	const Real JointState::getJointTorque() const
	{
		return _tau;
	}

	const SE3 & JointState::getJointT() const
	{
		return _T;
	}

	const SE3 & JointState::getJointAccumulatedT() const
	{
		return _accumulatedT;
	}

	const se3 & JointState::getJointScrew() const
	{
		return _screw;
	}

	const se3 & JointState::getJointScrewDot() const
	{
		return _screwDot;
	}

	const dse3 & JointState::getJointConstraintF() const
	{
		return _constraintF;
	}

	const Real JointState::getJointPos() const
	{
		return _q;
	}

	const Real JointState::getJointVel() const
	{
		return _qdot;
	}

	const Real JointState::getJointAcc() const
	{
		return _qddot;
	}

	bool JointState::checkJointInfoUpToDate(int infoIdx)
	{
		return infoIdx == (infoIdx & _jointInfoUpToDate);
	}

	void JointState::updateJointInfoUpToDate(int infoIdx, bool upToDate)
	{
		if (upToDate)
			_jointInfoUpToDate |= infoIdx;
		else
			_jointInfoUpToDate &= ~infoIdx;
	}

	void JointState::addJointPos(const Real q)
	{
		_q += q;
		updateJointInfoUpToDate(EXPOENTIAL, false);
	}

	void JointState::addJointVel(const Real qdot)
	{
		_qdot += qdot;
	}

	void JointState::addJointAcc(const Real qddot)
	{
		_qddot += qddot;
	}

	JointStatePtr JointState::copy() const
	{
		JointStatePtr clone(new JointState());
		clone->setJointTorque(_tau);
		clone->setJointT(_T);
		clone->setJointAccumulatedT(_accumulatedT);
		clone->setJointScrew(_screw);
		clone->setJointScrewDot(_screwDot);
		clone->setJointConstraintF(_constraintF);
		clone->setJointPos(_q);
		clone->setJointVel(_qdot);
		clone->setJointAcc(_qddot);
		return clone;
	}	

}