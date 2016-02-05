#pragma once

#include <vector>
#include <rovin/Dynamics/SerialOpenChain.h>

class efortRobot : public rovin::SerialOpenChain
{
public:
	efortRobot() : SerialOpenChain()
	{
		unsigned int DOF = 6;

		rovin::VectorX linkLength(6);
		linkLength << 0.504, 0.170, 0.780, 0.140, 0.760, 0.125;

		for (unsigned int i = 0; i < DOF + 1; i++)
		{
			addLink(rovin::LinkPtr(new rovin::Link()));
		}

		for (unsigned int i = 0; i < DOF; i++)
		{
			addMotorJoint(rovin::MotorJointPtr(new rovin::MotorJointPtr()));
		}
		addMate(0, rovin::SE3(rovin::Vector3(0, 0, linkLength(0))), rovin::SE3(rovin::Vector3(0, 0, 0)));
		addMate(1, rovin::SE3(rovin::SO3::RotX(-rovin::PI_HALF)*rovin::SO3::RotZ(-rovin::PI_HALF), rovin::Vector3(linkLength(1), 0, 0)), rovin::SE3(rovin::Vector3(0, 0, 0)));
		addMate(2, rovin::SE3(rovin::Vector3(linkLength(2), 0, 0)), rovin::SE3(rovin::Vector3(0, 0, 0)));
		addMate(3, rovin::SE3(rovin::SO3::RotX(-rovin::PI_HALF), rovin::Vector3(linkLength(3), linkLength(4), 0)), rovin::SE3(rovin::Vector3(0, 0, 0)));
		addMate(4, rovin::SE3(rovin::SO3::RotX(rovin::PI_HALF), rovin::Vector3(0, 0, 0)), rovin::SE3(rovin::Vector3(0, 0, 0)));
		addMate(5, rovin::SE3(rovin::SO3::RotX(-rovin::PI_HALF), rovin::Vector3(0, 0, 0)), rovin::SE3(rovin::Vector3(0, 0, 0)));
	}
};