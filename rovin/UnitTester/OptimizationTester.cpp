#include <rovin\Optimizer\NonlinearOptimization.h>
#include <rovin\Math\Function.h>
#include <rovin\EnergyOptimization\PTPOptimization.h>
#include <rovin\Renderer\OSG_simpleRender.h>

#include "efortRobot.h"

#include <iostream>
#include <fstream>
#include <conio.h>
#include <time.h>
#include <string.h>

using namespace rovin;
using namespace std;

class TestObjFunction;
class TestEqFunction;
class TestIneqFunction;

void saveMatrixX2txt(MatrixX in, std::string filename);
void saveVectorX2txt(VectorX in, std::string filename);
void calculateTorqueTrajectory(const MatrixX& q, const MatrixX& qdot, const MatrixX& qddot, MatrixX& torque);

std::string file = "C:/Users/crazy/Desktop/Time optimization/nloptMMA test/";

SerialOpenChainPtr robot(new efortRobot());

int main()
{
	unsigned int dof = robot->getNumOfJoint();

	StatePtr initState, finalState;
	initState = robot->makeState();
	finalState = robot->makeState();

	VectorX init_q(dof), init_qdot(dof), init_qddot(dof);
	VectorX final_q(dof), final_qdot(dof), final_qddot(dof);

	init_q << 1.0854, 1.02654, 0.798359, 2.97849, 1.50724, 1.45496;
	init_qdot.setZero(); init_qddot.setZero();
	final_q << -1.31465, 0.128787, -0.546992, 2.83671, -1.8583, 2.95029;
	final_qdot.setZero(); final_qddot.setZero();

	initState->setJointStatePos(init_q);initState->setJointStateVel(init_qdot);initState->setJointStateAcc(init_qddot);
	finalState->setJointStatePos(final_q);finalState->setJointStateVel(final_qdot);finalState->setJointStateAcc(final_qddot);

	///////////////////////////////// OPTIMIZATION /////////////////////////////////
	vector<bool> optJoint(robot->getNumOfJoint());
	optJoint[0] = optJoint[1] = optJoint[2] = true;
	Real tf = 2.0;
	int numOfOptCP = 4;
	int orderOfBSpline = 4;

	std::cout << "------- [NLOPT RESULT] -------" << endl;
	PTPOptimization PTPManagerNlopt(robot, optJoint, orderOfBSpline, numOfOptCP, 20, tf, initState, finalState, OptimizationType::nlopt);
	PTPManagerNlopt.generateTrajectory();

	std::cout << "------- [GCMMA RESULT] -------" << endl;
	PTPOptimization PTPManagerManualOpt(robot, optJoint, orderOfBSpline, numOfOptCP, 20, tf, initState, finalState, OptimizationType::GCMMA);
	PTPManagerManualOpt.generateTrajectory();

	///////////////////////////////// SAVE & RENDERING /////////////////////////////////
	bool renderingswi = false;

	int datanum = 2000;
	Real stepsize = (tf - 0.0) / 2000;
	Real t = 0.0;
	BSpline<-1, -1, -1> nloptSpline = PTPManagerNlopt._shared->_qSpline;
	BSpline<-1, -1, -1> MMASpline = PTPManagerManualOpt._shared->_qSpline;
	BSpline<-1, -1, -1> nloptqdotSpline = PTPManagerNlopt._shared->_qdotSpline;
	BSpline<-1, -1, -1> MMAqdotSpline = PTPManagerManualOpt._shared->_qdotSpline;
	BSpline<-1, -1, -1> nloptqddotSpline = PTPManagerNlopt._shared->_qddotSpline;
	BSpline<-1, -1, -1> MMAqddotSpline = PTPManagerManualOpt._shared->_qddotSpline;
	MatrixX nloptTraj(dof, datanum), MMATraj(dof, datanum);
	MatrixX nloptqdotTraj(dof, datanum), MMAqdotTraj(dof, datanum);
	MatrixX nloptqddotTraj(dof, datanum), MMAqddotTraj(dof, datanum);

	for (int i = 0; i < datanum; i++)
	{
		nloptTraj.col(i) = nloptSpline(t);
		MMATraj.col(i) = MMASpline(t);
		t += stepsize;
	}
	t = 0.0;
	for (int i = 0; i < datanum; i++)
	{
		nloptqdotTraj.col(i) = nloptqdotSpline(t);
		MMAqdotTraj.col(i) = MMAqdotSpline(t);
		t += stepsize;
	}
	t = 0.0;
	for (int i = 0; i < datanum; i++)
	{
		nloptqddotTraj.col(i) = nloptqddotSpline(t);
		MMAqddotTraj.col(i) = MMAqddotSpline(t);
		t += stepsize;
	}

	MatrixX nloptTorqueTraj(dof, datanum), MMATorqueTraj(dof, datanum);
	calculateTorqueTrajectory(nloptTraj, nloptqdotTraj, nloptqddotTraj, nloptTorqueTraj);
	calculateTorqueTrajectory(MMATraj, MMAqdotTraj, MMAqddotTraj, MMATorqueTraj);

	saveMatrixX2txt(nloptTraj, file + "nlopt/q.txt");
	saveMatrixX2txt(nloptqdotTraj, file + "nlopt/qdot.txt");
	saveMatrixX2txt(nloptqddotTraj, file + "nlopt/qddot.txt");
	saveMatrixX2txt(nloptTorqueTraj, file + "nlopt/torque.txt");

	saveMatrixX2txt(MMATraj, file + "MMA/q.txt");
	saveMatrixX2txt(MMAqdotTraj, file + "MMA/qdot.txt");
	saveMatrixX2txt(MMAqddotTraj, file + "MMA/qddot.txt");
	saveMatrixX2txt(MMATorqueTraj, file + "MMA/torque.txt");

	// constraint save
	VectorX qct(dof * 2), qdotct(dof * 2), qddotct(dof * 2), tct(dof * 2);
	for (int i = 0; i < dof; i++)
	{
		qct(i) = robot->getMotorJointPtr(i)->getLimitPosLower();
		qct(i + dof) = robot->getMotorJointPtr(i)->getLimitPosUpper();
		qdotct(i) = robot->getMotorJointPtr(i)->getLimitVelLower();
		qdotct(i + dof) = robot->getMotorJointPtr(i)->getLimitVelUpper();
		qddotct(i) = robot->getMotorJointPtr(i)->getLimitAccLower();
		qddotct(i + dof) = robot->getMotorJointPtr(i)->getLimitAccUpper();
		tct(i) = robot->getMotorJointPtr(i)->getLimitTorqueLower();
		tct(i + dof) = robot->getMotorJointPtr(i)->getLimitTorqueUpper();
	}

	saveVectorX2txt(qct, file + "qconstraint.txt");
	saveVectorX2txt(qdotct, file + "qdotconstraint.txt");
	saveVectorX2txt(qddotct, file + "qddotconstraint.txt");
	saveVectorX2txt(tct, file + "torqueconstraint.txt");

	if (renderingswi)
	{
		OSG_simpleRender renderer(*robot, *initState, 800, 800);
		shared_ptr<Points> sp, ep;
		shared_ptr<Line> l_nlopt, l_MMA;
		Real psize = 20.0;

		SE3 T(Vector3(0.0, 0.0, 0.12));
		initState->setJointStatePos(init_q); robot->solveForwardKinematics(*initState);
		Vector3 spos = (initState->getLinkStateSE3(dof)*T).getPosition();
		finalState->setJointStatePos(final_q); robot->solveForwardKinematics(*finalState);
		Vector3 epos = (finalState->getLinkStateSE3(dof)*T).getPosition();
		sp = shared_ptr<Points>(new Points);
		ep = shared_ptr<Points>(new Points);
		sp->push_back(osg::Vec3(spos(0), spos(1), spos(2)));
		ep->push_back(osg::Vec3(epos(0), epos(1), epos(2)));
		sp->setSize(psize); ep->setSize(psize);
		sp->setColor(0.0f, 0.0f, 1.0f); ep->setColor(1.0f, 0.0f, 0.0f);
		renderer.addGeometry(*sp);
		renderer.addGeometry(*ep);

		Vector3 nloptPos, MMAPos;
		l_nlopt = shared_ptr<Line>(new Line);
		l_MMA = shared_ptr<Line>(new Line);
		for (int i = 0; i < datanum; i++)
		{
			initState->setJointStatePos(nloptTraj.col(i));
			robot->solveForwardKinematics(*initState);
			nloptPos = (initState->getLinkStateSE3(dof)*T).getPosition();
			l_nlopt->push_back(osg::Vec3(nloptPos(0), nloptPos(1), nloptPos(2)));
		}
		for (int i = 0; i < datanum; i++)
		{
			initState->setJointStatePos(MMATraj.col(i));
			robot->solveForwardKinematics(*initState);
			MMAPos = (initState->getLinkStateSE3(dof)*T).getPosition();
			l_MMA->push_back(osg::Vec3(MMAPos(0), MMAPos(1), MMAPos(2)));
		}
		l_nlopt->setColor(0.0f, 0.0f, 1.0f);
		l_MMA->setColor(1.0f, 0.0f, 0.0f);
		renderer.addGeometry(*l_nlopt);
		renderer.addGeometry(*l_MMA);

		Real startT = 0.0, endT = tf;
		Real frame = 500;
		Real dT = 1.0 / frame;

		renderer.getViewer().realize();
		int count = 0;
		Real pastT = clock();

		while (1)
		{
			if ((clock() - pastT) / 1e+3 >= dT)
			{
				pastT = clock();

				if (count >= nloptTraj.cols())
				{
					count = 0;
				}
				initState->setJointStatePos(nloptTraj.col(count));
				robot->solveForwardKinematics(*initState);
				count += 5;
			}
			renderer.updateFrame();
		}
	}

	std::cout << "====== Program Complete ======" << endl;
	_getch();

	return 0;
}

void saveMatrixX2txt(MatrixX in, std::string filename)
{
	std::ofstream fout;
	fout.open(filename);

	for (int i = 0; i < in.cols(); i++)
	{
		for (int j = 0; j < in.rows(); j++)
			fout << in(j, i) << '\t';
		fout << std::endl;
	}

	fout.close();
}

void saveVectorX2txt(VectorX in, std::string filename)
{
	std::ofstream fout;
	fout.open(filename);

	for (int i = 0; i < in.size(); i++)
	{
		fout << in(i) << endl;
	}
	fout.close();
}

void calculateTorqueTrajectory(const MatrixX& q, const MatrixX& qdot, const MatrixX& qddot, MatrixX& torque)
{
	int dof = q.rows();
	int datanum = q.cols();
	StatePtr state = robot->makeState();

	for (int i = 0; i < datanum; i++)
	{
		state->setJointStatePos(q.col(i));
		state->setJointStateVel(qdot.col(i));
		state->setJointStateAcc(qddot.col(i));

		robot->solveInverseDynamics(*state);
		torque.col(i) = state->getJointStateTorque();
	}
}

class TestObjFunction : public Function
{
public:
	TestObjFunction() {}

	VectorX func(const VectorX& x) const
	{
		VectorX val(1);
		val(0) = x(0)*x(0) + x(1)*x(1);
		return val;
	}
};

class TestEqFunction : public Function
{
public:
	TestEqFunction() {}

	VectorX func(const VectorX& x) const
	{
		VectorX val(1);
		val(0) = x(0) + x(1) - 5;
		return val;
	}
};

class TestIneqFunction : public Function
{
public:
	TestIneqFunction() {}

	VectorX func(const VectorX& x) const
	{
		VectorX val(1);
		val(0) = -x(1) + 3;
		return val;
	}
};