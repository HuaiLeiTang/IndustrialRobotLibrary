#include "MotorJoint.h"

namespace rovin {
	MotorJoint::MotorJoint() : _rotorinertia(0), _resistance(0), _inductance(0), _gearRatio(0),
		_springConstant(0), _damperConstant(0), _ViscousFrictionConstant(0), _CoulombFrictionConstant(0)
	{
		_axis.setZero();
		// z-axis (w,v) = (0, 0, 1, 0, 0, 0)
		_axis[2] = 1; 

		_LimitPosLower = std::numeric_limits<Real>::min();
		_LimitPosUpper = std::numeric_limits<Real>::max();

		_LimitVelLower = std::numeric_limits<Real>::min();
		_LimitVelUpper = std::numeric_limits<Real>::max();

		_LimitAccLower = std::numeric_limits<Real>::min();
		_LimitAccUpper = std::numeric_limits<Real>::max();

		_LimitJerkLower = std::numeric_limits<Real>::min();
		_LimitJerkUpper = std::numeric_limits<Real>::max();

		_LimitTorqueLower = std::numeric_limits<Real>::min();
		_LimitTorqueUpper = std::numeric_limits<Real>::max();
	}

	MotorJoint::MotorJoint(const Vector6 & axis) : _axis(axis), _rotorinertia(0), _resistance(0), _inductance(0), _gearRatio(0),
		_springConstant(0), _damperConstant(0), _ViscousFrictionConstant(0), _CoulombFrictionConstant(0)
	{
		_LimitPosLower = std::numeric_limits<Real>::min();
		_LimitPosUpper = std::numeric_limits<Real>::max();

		_LimitVelLower = std::numeric_limits<Real>::min();
		_LimitVelUpper = std::numeric_limits<Real>::max();

		_LimitAccLower = std::numeric_limits<Real>::min();
		_LimitAccUpper = std::numeric_limits<Real>::max();

		_LimitJerkLower = std::numeric_limits<Real>::min();
		_LimitJerkUpper = std::numeric_limits<Real>::max();

		_LimitTorqueLower = std::numeric_limits<Real>::min();
		_LimitTorqueUpper = std::numeric_limits<Real>::max();
	}

	MotorJoint::~MotorJoint() {}

	void MotorJoint::setRotorInertia(const Real RI)
	{
		_rotorinertia = RI;
	}

	void MotorJoint::setResistance(const Real R)
	{
		_resistance = R;
	}

	void MotorJoint::setInductance(const Real In)
	{
		_inductance = In;
	}

	void MotorJoint::setGearRatio(const Real Gear)
	{
		_gearRatio = Gear;
	}

	void MotorJoint::setSpringConstant(const Real Spring)
	{
		_springConstant = Spring;
	}

	void MotorJoint::setDamperConstant(const Real Damper)
	{
		_damperConstant = Damper;
	}

	void MotorJoint::setViscousFrictionConstant(const Real VFriction)
	{
		_ViscousFrictionConstant = VFriction;
	}

	void MotorJoint::setCoulombFrictionConstant(const Real CFriction)
	{
		_CoulombFrictionConstant = CFriction;
	}

	bool MotorJoint::setLimitPos(const Real lower, const Real upper)
	{
		_LimitPosLower = lower;
		_LimitPosUpper = upper;
		return true;
	}

	bool MotorJoint::setLimitVel(const Real lower, const Real upper)
	{
		_LimitVelLower = lower;
		_LimitVelUpper = upper;
		return true;
	}

	bool MotorJoint::setLimitAcc(const Real lower, const Real upper)
	{
		_LimitAccLower = lower;
		_LimitAccUpper = upper;
		return true;
	}

	bool MotorJoint::setLimitJerk(const Real lower, const Real upper)
	{
		_LimitJerkLower = lower;
		_LimitJerkUpper = upper;
		return true;
	}

	bool MotorJoint::setLimitTorque(const Real lower, const Real upper)
	{
		_LimitTorqueLower = lower;
		_LimitTorqueUpper = upper;
		return true;
	}

	void MotorJoint::setAxis(const Vector3 & axis)
	{
		_axis[0] = axis[0];
		_axis[1] = axis[1];
		_axis[2] = axis[2];
	}

	void MotorJoint::setAxis(const Vector6 & axis)
	{
		_axis = axis;
	}

	const Real MotorJoint::getRotorInertia() const
	{
		return _rotorinertia;
	}

	const Real MotorJoint::getResistance() const
	{
		return _resistance;
	}

	const Real MotorJoint::getInductance() const
	{
		return _inductance;
	}

	const Real MotorJoint::getGearRatio() const
	{
		return _gearRatio;
	}

	const Real MotorJoint::getSpringConstant() const
	{
		return _springConstant;
	}

	const Real MotorJoint::getDamperConstant() const
	{
		return _damperConstant;
	}

	const Real MotorJoint::getViscousFrictionConstant() const
	{
		return _ViscousFrictionConstant;
	}

	const Real MotorJoint::getCoulombFrictionConstant() const
	{
		return _CoulombFrictionConstant;
	}

	const Real MotorJoint::getLimitPosLower() const
	{
		return _LimitPosLower;
	}

	const Real MotorJoint::getLimitPosUpper() const
	{
		return _LimitPosUpper;
	}

	const Real MotorJoint::getLimitVelLower() const
	{
		return _LimitVelLower;
	}

	const Real MotorJoint::getLimitVelUpper() const
	{
		return _LimitVelUpper;
	}

	const Real MotorJoint::getLimitAccLower() const
	{
		return _LimitAccLower;
	}

	const Real MotorJoint::getLimitAccUpper() const
	{
		return _LimitAccUpper;
	}

	const Real MotorJoint::getLimitJerkLower() const
	{
		return _LimitJerkLower;
	}

	const Real MotorJoint::getLimitJerkUpper() const
	{
		return _LimitJerkUpper;
	}

	const Real MotorJoint::getLimitTorqueLower() const
	{
		return _LimitTorqueLower;
	}

	const Real MotorJoint::getLimitTorqueUpper() const
	{
		return _LimitTorqueUpper;
	}

	const Vector6 & MotorJoint::getAxis() const
	{
		return _axis;
	}

	MotorJointPtr MotorJoint::copy() const
	{
		MotorJointPtr clone(new MotorJoint(_axis));
		clone->setRotorInertia(_rotorinertia);
		clone->setResistance(_resistance);
		clone->setInductance(_inductance);
		clone->setGearRatio(_gearRatio);
		clone->setSpringConstant(_springConstant);
		clone->setDamperConstant(_damperConstant);
		clone->setViscousFrictionConstant(_ViscousFrictionConstant);
		clone->setCoulombFrictionConstant(_CoulombFrictionConstant);
		clone->setLimitPos(_LimitPosLower, _LimitPosUpper);
		clone->setLimitVel(_LimitVelLower, _LimitVelUpper);
		clone->setLimitAcc(_LimitAccLower, _LimitAccUpper);
		clone->setLimitJerk(_LimitJerkLower, _LimitJerkUpper);
		clone->setLimitTorque(_LimitTorqueLower, _LimitTorqueUpper);
		return clone;
	}

}