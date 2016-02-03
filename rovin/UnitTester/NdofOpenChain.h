#pragma once

#include <vector>
#include <rovin/Dynamics/SerialOpenChain.h>

template <unsigned int dof>
class NdofOpenChain : public rovin::SerialOpenChain
{
public:
	NdofOpenChain() : SerialOpenChain()
	{
		for (unsigned int i = 0; i < dof + 1; i++)
		{
			rovin::LinkPtr link_tmp = rovin::LinkPtr(new rovin::Link());
			link_tmp->addDrawingGeomtryInfo(std::shared_ptr< rovin::Box >(new Box(3, 3, 15)));

			addLink(link_tmp);
		}

		for (unsigned int i = 0; i < dof; i++)
		{
			addMotorJoint(rovin::MotorJointPtr(new rovin::MotorJoint()));
			addMate(i,
				rovin::SE3(rovin::SO3::RotY(-rovin::PI_HALF), rovin::Vector3(10, 0, 0)),
				rovin::SE3(rovin::SO3::RotY(-rovin::PI_HALF), rovin::Vector3(-10, 0, 0)));
		}

		completeAssembling();
	}
};