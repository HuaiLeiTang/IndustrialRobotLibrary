#include <cstdio>
#include <conio.h>
#include <ctime>`
#include <cstdlib>
#include <iostream>
#include <rovin/Renderer/OSG_simpleRender.h>

#include "NdofOpenChain.h"
#include <rovin/Utils/Diagnostic.h>

using namespace std;
using namespace rovin;

const unsigned int DOF = 3;

int main()
{
	srand(time(NULL));

	NdofOpenChain< DOF > robot;
	StatePtr state = robot.makeState();
	SE3 GoalT;

	VectorX q(DOF), qdot(DOF), qddot(DOF);

	q << 3, 2, 1;
	qdot << 3, 2, 1;
	qddot << 3, 2, 1;

	state->setJointStatePos(q);
	state->setJointStateVel(qdot);
	state->setJointStateAcc(qddot);

	//state->setJointStatePos(0, 1);
	//state->setJointStatePos(1, 2);
	//state->setJointStatePos(2, PI_HALF);

	//state->setJointStateVel(0, 1);
	//state->setJointStateVel(1, 2);
	//state->setJointStateVel(2, 3);

	//state->setJointStateAcc(0, 1);
	//state->setJointStateAcc(1, 2);
	//state->setJointStateAcc(2, 3);

	PERFORM_TEST(
		q.setRandom();
	qdot.setRandom();
	qddot.setRandom();

	state->setJointStatePos(q);
	state->setJointStateVel(qdot);
	state->setJointStateAcc(qddot);

	robot.solveInverseDynamics(*state);,
		1e+6
		);

	cout << "q = " << state->getJointStatePos() << endl;
	cout << "qdot = " << state->getJointStateVel() << endl;
	cout << "qddot = " << state->getJointStateAcc() << endl;

	//robot.solveInverseDynamics(*state);

	cout << "VEL = " << endl << state->getLinkStateVel(3) << endl;
	cout << "ACC = " << endl << state->getLinkStateAcc(3) << endl;

	cout << "tau(0) = " << state->getJointStateTorque(0) << endl;
	cout << "tau(1) = " << state->getJointStateTorque(1) << endl;
	cout << "tau(2) = " << state->getJointStateTorque(2) << endl;
	cout << endl;

	//cout << "q(0) = " << state->getJointStatePos(0) << endl;
	//cout << "q(1) = " << state->getJointStatePos(1) << endl;
	//cout << "q(2) = " << state->getJointStatePos(2) << endl;
	//cout << endl;

	//robot.solveForwardKinematics(*state);

	//cout << "T(3) = " << endl << (GoalT = state->getLinkStateSE3(3)) << endl;
	//cout << endl;

	//state->setJointStatePos(0, 0);
	//state->setJointStatePos(1, 0);
	//state->setJointStatePos(2, 0);

	//robot.solveInverseKinematics(*state, GoalT);

	//cout << "q(0) = " << state->getJointStatePos(0) << endl;
	//cout << "q(1) = " << state->getJointStatePos(1) << endl;
	//cout << "q(2) = " << state->getJointStatePos(2) << endl;
	//cout << endl;

	//robot.solveForwardKinematics(*state);

	//cout << "T(3) = " << endl << (GoalT = state->getLinkStateSE3(3)) << endl;
	//cout << endl;

	_getch();

	return 0;
}