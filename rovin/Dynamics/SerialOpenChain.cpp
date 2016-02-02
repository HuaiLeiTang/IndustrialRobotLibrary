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
		StatePtr state(new State(_motorJointPtr.size()));
		return state;
	}

	/*!
	* \brief Serial Open Chain kinematics functions
	*/
	void SerialOpenChain::solveForwardKinematics(State & state, JOINT_KINEMATICS_OPTION option)
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