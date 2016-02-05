/*!
 *	\file	MotorJoint.h
 *	\date	2016.01.22
 *	\author	Youngsuk (crazyhys@gmail.com)
 *	\brief	MotorJoin class
  *         description of the industrial robot motors
*/

#pragma once

#include <string>
#include <memory>

#include <Eigen\Dense>
#include <rovin\Math\LieGroup.h>
#include <rovin\Math\Constant.h>

namespace rovin{

	class MotorJoint;
	typedef std::shared_ptr<MotorJoint> MotorJointPtr;

	class MotorJoint
	{

	private:
		/*!
		* \brief MotorJoint class member variable
		*/

		// Motor Properties
		Real _rotorinertia;
		Real _resistance;
		Real _inductance;
		Real _gearRatio;
		Real _motorConstant;
		Real _backEMFConstant;

		Real _springConstant;
		Real _damperConstant;
		Real _ViscousFrictionConstant;
		Real _CoulombFrictionConstant;
		
		// Motor Sensor
		//class Encoder; // »ý°¢Á»..

		// Motor limit Values	
		Real _LimitPosLower;
		Real _LimitPosUpper;
		Real _LimitVelLower;
		Real _LimitVelUpper;
		Real _LimitAccLower;
		Real _LimitAccUpper;
		Real _LimitJerkLower;
		Real _LimitJerkUpper;
		Real _LimitTorqueLower;
		Real _LimitTorqueUpper;

		// Joint member values
		Vector6 _axis;

	public:
		/*!
		* \brief MotorJoint class member functions
		*/

		// constructor & destructor
		MotorJoint();
		MotorJoint(const Vector6& axis);
		~MotorJoint();

		// set-function
		void setRotorInertia(const Real RI);
		void setResistance(const Real R);
		void setInductance(const Real In);
		void setGearRatio(const Real Gear);
		void setMotorConstant(const Real MC);
		void setBackEMFConstant(const Real BEMFC);
		void setSpringConstant(const Real Spring);
		void setDamperConstant(const Real Damper);
		void setViscousFrictionConstant(const Real VFriction);
		void setCoulombFrictionConstant(const Real CFriction);

		bool setLimitPos(const Real lower, const Real upper);
		bool setLimitVel(const Real lower, const Real upper);
		bool setLimitAcc(const Real lower, const Real upper);
		bool setLimitJerk(const Real lower, const Real upper);
		bool setLimitTorque(const Real lower, const Real upper);

		void setAxis(const Vector3& axis); // when only setting in angular velocity
		void setAxis(const Vector6& axis);

		// get-function
		const Real getRotorInertia() const;
		const Real getResistance() const;
		const Real getInductance() const;
		const Real getGearRatio() const;
		const Real getMotorConstant() const;
		const Real getBackEMFConstant() const;
		const Real getSpringConstant() const;
		const Real getDamperConstant() const;
		const Real getViscousFrictionConstant() const;
		const Real getCoulombFrictionConstant() const;

		const Real getLimitPosLower() const;
		const Real getLimitPosUpper() const;
		const Real getLimitVelLower() const;
		const Real getLimitVelUpper() const;
		const Real getLimitAccLower() const;
		const Real getLimitAccUpper() const;
		const Real getLimitJerkLower() const;
		const Real getLimitJerkUpper() const;
		const Real getLimitTorqueLower() const;
		const Real getLimitTorqueUpper() const;

		const Vector6& getAxis() const;

		// deep-copy
		MotorJointPtr copy() const;

	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW;


	};

}

