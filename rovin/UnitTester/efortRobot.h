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
			addMotorJoint(rovin::MotorJointPtr(new rovin::MotorJoint()));
		}

		addMate(0, rovin::SE3(rovin::Vector3(0, 0, linkLength(0))).inverse(), rovin::SE3(rovin::Vector3(0, 0, 0)));
		addMate(1, rovin::SE3(rovin::SO3::RotX(-rovin::PI_HALF)*rovin::SO3::RotZ(-rovin::PI_HALF), rovin::Vector3(linkLength(1), 0, 0)).inverse(), rovin::SE3(rovin::Vector3(0, 0, 0)));
		addMate(2, rovin::SE3(rovin::Vector3(linkLength(2), 0, 0)).inverse(), rovin::SE3(rovin::Vector3(0, 0, 0)));
		addMate(3, rovin::SE3(rovin::SO3::RotX(-rovin::PI_HALF), rovin::Vector3(linkLength(3), linkLength(4), 0)).inverse(), rovin::SE3(rovin::Vector3(0, 0, 0)));
		addMate(4, rovin::SE3(rovin::SO3::RotX(rovin::PI_HALF), rovin::Vector3(0, 0, 0)).inverse(), rovin::SE3(rovin::Vector3(0, 0, 0)));
		addMate(5, rovin::SE3(rovin::SO3::RotX(-rovin::PI_HALF), rovin::Vector3(0, 0, 0)).inverse(), rovin::SE3(rovin::Vector3(0, 0, 0)));


		rovin::MatrixX parameter(10, 6);
		parameter << 1.0000, -0.2726, 0.4164, 0.0379, 0.0379, 0.0379
			, 1.0000, 72.0862, 6.3893, 0.3574, -0.1061, -0.0759
			, 1.0000, -0.0000, 23.4120, 0.0124, 0.1181, 0.0256
			, 1.0000, 1.0000, -0.4166, -1.0914, -0.0124, 0.1181
			, 1.0000, 14.9404, 0.1093, 1.6195, 0.8195, 0.8034
			, 1.0000, -14.8427, -0.0116, 0.1966, 0.4121, 0.8312
			, 0.0478, 9.4016, 1.7068, -0.3880, 0.8152, 0.5500
			, 1.0000, -1.2929, 5.4276, -0.1850, 0.1343, -0.0108
			, 1.0000, 0.0320, 0.6453, 0.1541, -0.1275, 0.0534
			, 1.0000, -4.4974, -0.1622, -0.1995, -0.0834, 0.0317;

		rovin::VectorX m = parameter.row(0);
		rovin::MatrixX mx = parameter.block(1, 0, 3, 6);
		rovin::MatrixX inr_info = parameter.block(4, 0, 6, 6).transpose();
		rovin::Matrix3 I;
		std::vector<rovin::Inertia, Eigen::aligned_allocator<rovin::Inertia>> G(6);
		rovin::Vector3 p;
		for (int i = 0; i < 6; i++)
		{
			I(0, 0) = inr_info(i, 0);
			I(1, 1) = inr_info(i, 1);
			I(2, 2) = inr_info(i, 2);
			I(0, 1) = I(1, 0) = inr_info(i, 3);
			I(0, 2) = I(2, 0) = inr_info(i, 4);
			I(1, 2) = I(2, 1) = inr_info(i, 5);

			//std::cout << I << std::endl;
			//std::cout << p << std::endl;
			p = -mx.col(i);
			G[i] = rovin::Inertia(I, p, m(i));
		}

		for (unsigned int i = 0; i < DOF; i++)
		{
			getLinkPtr(i + 1)->setInertia(G[i]);
		}

		rovin::VectorX L(6);
		rovin::VectorX R(6);
		rovin::VectorX kt(6);
		rovin::VectorX kb(6);
		rovin::VectorX J(6);
		rovin::VectorX gearRatio(6);
		L << 3.5, 3.5, 5.2, 8, 8, 8;
		L *= 1e-3;

		R << 0.58, 0.58, 0.8, 2.9, 2.9, 7.5;
		kt << 0.73, 0.73, 0.5, 0.4, 0.4, 0.39;
		kb = kt;
		J << 1.06, 1.06, 0.13, 0.044, 0.044, 0.027;
		J *= 1e-3;
		gearRatio << 147, 153, 153, 76.95, 80, 51;
		for (int i = 0; i < 6; i++)
		{
			getMotorJointPtr(i)->setInductance(L(i));
			getMotorJointPtr(i)->setResistance(R(i));
			getMotorJointPtr(i)->setMotorConstant(kt(i));
			getMotorJointPtr(i)->setBackEMFConstant(kb(i));
			getMotorJointPtr(i)->setRotorInertia(J(i));
			getMotorJointPtr(i)->setGearRatio(gearRatio(i));
		}

		// set constraints

		rovin::VectorX taumax(6);
		rovin::VectorX taumin(6);
		rovin::VectorX qmax(6);
		rovin::VectorX qmin(6);
		rovin::VectorX qdotmax(6);
		rovin::VectorX qdotmin(6);
		rovin::VectorX qddotmax(6);
		rovin::VectorX qddotmin(6);
		rovin::VectorX qdddotmax(6);
		rovin::VectorX qdddotmin(6);
		rovin::VectorX kv(6);
		rovin::VectorX kc(6);
		rovin::VectorX smax(6);
		rovin::VectorX imax(6);

		qmax << 175, 90, 70, 180, 135, 360;
		qmax *= rovin::PI / 180;

		qmin << 175, 100, 145, 180, 135, 360;
		qmin *= -rovin::PI / 180;
		qdotmax << 100, 80, 140, 290, 290, 440;
		qdotmax *= rovin::PI / 180;

		qddotmax = 5 * qdotmax;             // user defined
		qddotmin = -qddotmax;
		qdddotmax = 300 * rovin::VectorX::Ones(6);          // user defined
		qdddotmin = -qdddotmax;
		taumax = 3000 * rovin::VectorX::Ones(6);          // user defined

		kv << 105.4116, 107.5105, 8.1031, 6.5430, 11.3551, 4.2050;
		kc << 92.3012, 117.3703, 57.8389, 10.0524, 24.6819, 20.32;
		smax << 3000, 3000, 4500, 4500, 4500, 4500;
		smax *= 2 * rovin::PI / 60;    // smax = [5000, 5000, 5000, 5000, 5000, 50000] * 2 * pi / 60,    % motor maximum
		imax << 39, 39, 20, 12, 12, 7.2;

		rovin::VectorX temptau = (imax.cwiseProduct(gearRatio)).cwiseProduct(kt);
		rovin::VectorX tempvel = smax.cwiseQuotient(gearRatio);
		taumax = taumax.cwiseMin(temptau); // correct
		qdotmax = qdotmax.cwiseMin(tempvel); // correct
		//taumax = taumax.cwiseMax(temptau);
		//qdotmax = qdotmax.cwiseMax(tempvel);
		qdotmin = -qdotmax;
		taumin = -taumax;
		for (unsigned int i = 0; i < 6; i++)
		{
			//getMotorJointPtr(i)->setDamperConstant(kv(i));
			//getMotorJointPtr(i)->setCoulombFrictionConstant(kc(i));
			getMotorJointPtr(i)->setLimitPos(qmin(i), qmax(i));
			getMotorJointPtr(i)->setLimitVel(qdotmin(i), qdotmax(i));
			getMotorJointPtr(i)->setLimitAcc(qddotmin(i), qddotmax(i));
			getMotorJointPtr(i)->setLimitJerk(qdddotmin(i), qdddotmax(i));
			getMotorJointPtr(i)->setLimitTorque(taumin(i), taumax(i));
		}

		completeAssembling();

		rovin::StatePtr efortState = makeState();
		rovin::VectorX q(6);
		q.setZero();
		efortState->setJointStatePos(q);
		solveForwardKinematics(*efortState);
		rovin::Vector4	orange(254 / 255.0, 193 / 255.0, 27 / 255.0, 1.0),
			black(55 / 255.0, 55 / 255.0, 55 / 255.0, 1.0),
			white(200 / 255.0, 200 / 255.0, 200 / 255.0, 1.0);
		//for (unsigned int i = 0; i < _links.size(); i++)
		//{
		//	shared_ptr<Box> boxShape(new Box(0.2, 0.2, 0.2));
		//	this->getLinkPtr(i)->addDrawingShapes(boxShape);
		//}
		for (unsigned int i = 0; i < 1; i++)
		{
			std::shared_ptr<rovin::Mesh> STL_file(new rovin::Mesh(std::string("../Data/CAD/efort_robot/LINK0_0") + std::to_string(i + 1) + std::string(".STL")));
			STL_file->setFrame(efortState->getLinkStateSE3(0).inverse());
			STL_file->setDimension(1);
			STL_file->setColor(orange);
			getLinkPtr(0)->addDrawingGeomtryInfo(STL_file);
		}
		for (unsigned int i = 0; i < 6; i++)
		{
			std::shared_ptr<rovin::Mesh> STL_file(new rovin::Mesh(std::string("../Data/CAD/efort_robot/LINK1_0") + std::to_string(i + 1) + std::string(".STL")));
			STL_file->setFrame(efortState->getLinkStateSE3(1).inverse());
			STL_file->setDimension(1);
			if (i == 1 || i == 2 || i == 3 || i == 4)
				STL_file->setColor(black);
			else
				STL_file->setColor(orange);
			getLinkPtr(1)->addDrawingGeomtryInfo(STL_file);
		}
		for (unsigned int i = 0; i < 1; i++)
		{
			std::shared_ptr<rovin::Mesh> STL_file(new rovin::Mesh(std::string("../Data/CAD/efort_robot/LINK2_0") + std::to_string(i + 1) + std::string(".STL")));
			STL_file->setFrame(efortState->getLinkStateSE3(2).inverse());
			STL_file->setDimension(1);
			STL_file->setColor(orange);
			getLinkPtr(2)->addDrawingGeomtryInfo(STL_file);
		}
		for (unsigned int i = 0; i < 7; i++)
		{
			std::shared_ptr<rovin::Mesh> STL_file(new rovin::Mesh(std::string("../Data/CAD/efort_robot/LINK3_0") + std::to_string(i + 1) + std::string(".STL")));
			STL_file->setFrame(efortState->getLinkStateSE3(3).inverse());
			STL_file->setDimension(1);
			if (i == 1 || i == 2 || i == 3 || i == 4 || i == 5)
				STL_file->setColor(black);
			else
				STL_file->setColor(orange);
			getLinkPtr(3)->addDrawingGeomtryInfo(STL_file);
		}
		for (unsigned int i = 0; i < 8; i++)
		{
			std::shared_ptr<rovin::Mesh> STL_file(new rovin::Mesh(std::string("../Data/CAD/efort_robot/LINK4_0") + std::to_string(i + 1) + std::string(".STL")));
			STL_file->setFrame(efortState->getLinkStateSE3(4).inverse());
			STL_file->setDimension(1);
			if (i == 0)
				STL_file->setColor(black);
			else
				STL_file->setColor(orange);
			getLinkPtr(4)->addDrawingGeomtryInfo(STL_file);
		}
		for (unsigned int i = 0; i < 3; i++)
		{
			std::shared_ptr<rovin::Mesh> STL_file(new rovin::Mesh(std::string("../Data/CAD/efort_robot/LINK5_0") + std::to_string(i + 1) + std::string(".STL")));
			STL_file->setFrame(efortState->getLinkStateSE3(5).inverse());
			STL_file->setDimension(1);
			STL_file->setColor(orange);
			getLinkPtr(5)->addDrawingGeomtryInfo(STL_file);
		}
		for (unsigned int i = 0; i < 1; i++)
		{
			std::shared_ptr<rovin::Mesh> STL_file(new rovin::Mesh(std::string("../Data/CAD/efort_robot/LINK6_0") + std::to_string(i + 1) + std::string(".STL")));
			STL_file->setFrame(efortState->getLinkStateSE3(6).inverse());
			STL_file->setDimension(1);
			if (i == 0)
				STL_file->setColor(black);
			else
				STL_file->setColor(orange);
			getLinkPtr(6)->addDrawingGeomtryInfo(STL_file);
		}

		//sh edit
		std::shared_ptr<rovin::Box> Box_Assemble(new rovin::Box());
		Box_Assemble->setFrame(rovin::SE3(rovin::Vector3(0, 0, linkLength(5) + 0.04)));
		Box_Assemble->setDimension(0.06, 0.06, 0.08);
		Box_Assemble->setColor(white);
		getLinkPtr(6)->addDrawingGeomtryInfo(Box_Assemble);

		std::shared_ptr<rovin::Cylinder> Cylinder_Assemble(new rovin::Cylinder());
		Cylinder_Assemble->setFrame(rovin::SE3(rovin::Vector3(0.0, 0.0, linkLength(5))) * rovin::SE3(rovin::SO3::RotZ(rovin::PI)*rovin::SO3::RotY(-rovin::PI_HALF), rovin::Vector3(0.04, 0.0, 0.049)));
		Cylinder_Assemble->setDimension(0.015, 0.12);
		Cylinder_Assemble->setColor(white);
		getLinkPtr(6)->addDrawingGeomtryInfo(Cylinder_Assemble);

		std::shared_ptr<rovin::Cylinder> Cylinder_Assemble2(new rovin::Cylinder());
		Cylinder_Assemble2->setFrame(rovin::SE3(rovin::Vector3(0.0, 0.0, linkLength(5))) * rovin::SE3(rovin::SO3::RotZ(rovin::PI)*rovin::SO3::RotY(-rovin::PI_HALF), rovin::Vector3(0.1, 0.0, 0.049)));
		Cylinder_Assemble2->setDimension(0.002, 0.0525 * 2);
		Cylinder_Assemble2->setColor(white);
		getLinkPtr(6)->addDrawingGeomtryInfo(Cylinder_Assemble2);

		// wy edit
		// end effector: welding gun
		//std::shared_ptr<rovin::Mesh> STL_file_WD(new rovin::Mesh(std::string("../Data/CAD/efort_robot/welding_gun_82W.STL")));
		//rovin::SE3 TlastLinkToEndeffector = rovin::SE3(rovin::Vector3(0.0, 0.0, 0.125)) * rovin::SE3(rovin::SO3::RotZ(rovin::PI)*rovin::SO3::RotY(-rovin::PI_HALF), rovin::Vector3(0.1525, 0.0, 0.0490));
		//STL_file_WD->setFrame(TlastLinkToEndeffector.inverse() * efortState->getLinkStateSE3(6).inverse());
		//STL_file_WD->setDimension(0.001);
		//STL_file_WD->setColor(white);
		//getLinkPtr(6)->addDrawingGeomtryInfo(STL_file_WD);

		//for (unsigned int i = 0; i < 3; i++)
		//{
		//	std::shared_ptr<rovin::Mesh> STL_file(new rovin::Mesh(std::string("../Data/CAD/efort_robot/tooltip") + std::to_string(i + 1) + std::string(".STL")));
		//	STL_file->setFrame(efortState->getLinkStateSE3(6).inverse() * rovin::SE3(rovin::Vector3(1.104, 0.0, 1.5765))*TlastLinkToEndeffector.inverse());
		//	STL_file->setDimension(0.001);
		//	STL_file->setColor(white);
		//	getLinkPtr(6)->addDrawingGeomtryInfo(STL_file);
		//}
	}
};