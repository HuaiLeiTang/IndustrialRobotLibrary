#include <rovin\Optimizer\NonlinearOptimization.h>
#include <rovin\Math\Function.h>
#include <rovin\EnergyOptimization\PTPOptimization.h>
#include <rovin\EnergyOptimization\PTPWayPointOptimization.h>

#include "efortRobot.h"

#include <iostream>
#include <conio.h>

#include <string>
#include <memory>

using namespace rovin;
using namespace std;

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

class human
{
public:
	string name;
public:
	human() { name = "NONE"; }
	human(string _name) : name(_name) {}
	void getname()
	{
		cout << "Human : " << name << endl;
	}
};

class Hong : public human
{
public:
	string name;
	int age;
public:
	Hong() { name = "HONG NONE"; age = 0; }
	Hong(string _name, int _age) : name(_name), age(_age) {}
	void getname()
	{
		cout << "Hong : " << name << endl;
	}
	void getage()
	{
		cout << "Hong : " << age << endl;
	}
};

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

	initState->setJointStatePos(init_q);
	initState->setJointStateVel(init_qdot);
	initState->setJointStateAcc(init_qddot);

	finalState->setJointStatePos(final_q);
	finalState->setJointStateVel(final_qdot);
	finalState->setJointStateAcc(final_qddot);

	vector<bool> optJoint(robot->getNumOfJoint());
	optJoint[0] = optJoint[1] = optJoint[2] = true;

	int orderOfspline = 4;
	int numOfCP = 4;
	//int numOfCP = 6;
	int numOfGQSample = 20;
	Real tf = 2.0;
	
	PTPOptimization PTPManager(robot, optJoint, orderOfspline, numOfCP, numOfGQSample, tf, initState, finalState);
	PTPManager.generateTrajectory();

	PTPWayPointOptimization PTPWayPointManager(robot, optJoint, orderOfspline, numOfCP, numOfGQSample, tf, initState, finalState, PTPWayPointOptimization::INITIAL_ALL_FINAL_POS);
	PTPWayPointManager.generateTrajectory();

	VectorX CP = PTPWayPointManager._shared->_qSpline.getControlPoints().row(0);
	VectorX knot = PTPWayPointManager._knot;

	cout << "CP" << endl;
	cout << CP << endl;
	cout << "knot" << endl;
	cout << knot << endl;


	BSpline<-1, -1, -1> tmpSpline(knot, CP.transpose());
	BSpline<-1, -1, -1> tmpSplinedot = tmpSpline.derivative();
	BSpline<-1, -1, -1> tmpSplineddot = tmpSplinedot.derivative();
	
	cout << tmpSpline.fval(tf- 0.00000001) << endl;
	cout << tmpSplinedot.fval(tf - 0.00000001) << endl;
	cout << tmpSplineddot.fval(tf - 0.00000001) << endl;



	_getch();

	return 0;
}