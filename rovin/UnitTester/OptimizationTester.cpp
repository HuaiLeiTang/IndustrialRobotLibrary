#include <rovin\Optimizer\NonlinearOptimization.h>
#include <rovin\Math\Function.h>
#include <rovin\EnergyOptimization\PTPOptimization.h>
#include <rovin\Renderer\OSG_simpleRender.h>

#include "efortRobot.h"

#include <iostream>
#include <conio.h>
#include <time.h>

using namespace rovin;
using namespace std;

class TestObjFunction;
class TestEqFunction;
class TestIneqFunction;

int main()
{
	SerialOpenChainPtr robot(new efortRobot());
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
	int numOfOptCP = 6;
	int orderOfBSpline = 4;

	cout << "------- [NLOPT RESULT] -------" << endl;
	PTPOptimization PTPManagerNlopt(robot, optJoint, orderOfBSpline, numOfOptCP, 20, tf, initState, finalState, OptimizationType::nlopt);
	PTPManagerNlopt.generateTrajectory();

	cout << "------- [GCMMA RESULT] -------" << endl;
	PTPOptimization PTPManagerManualOpt(robot, optJoint, orderOfBSpline, numOfOptCP, 20, tf, initState, finalState, OptimizationType::GCMMA);
	PTPManagerManualOpt.generateTrajectory();

	///////////////////////////////// RENDERING /////////////////////////////////
	bool renderingswi = false;

	BSpline<-1, -1, -1> nloptSpline = PTPManagerNlopt._shared->_qSpline;
	BSpline<-1, -1, -1> MMASpline = PTPManagerManualOpt._shared->_qSpline;
	int datanum = 2000;
	Real stepsize = (tf - 0.0) / 2000;
	Real t = 0.0;
	MatrixX nloptTraj(dof, datanum), MMATraj(dof, datanum);
	for (int i = 0; i < datanum; i++)
	{
		nloptTraj.col(i) = nloptSpline(t);
		MMATraj.col(i) = MMASpline(t);
		t += stepsize;
	}

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

	_getch();

	return 0;
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