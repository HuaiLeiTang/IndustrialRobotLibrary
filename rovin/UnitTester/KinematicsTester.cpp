#include <cstdio>
#include <conio.h>
#include <ctime>`
#include <cstdlib>
#include <rovin/Renderer/OSG_simpleRender.h>

#include "NdofOpenChain.h"
#include <rovin/Utils/Diagnostic.h>

using namespace std;
using namespace rovin;

const unsigned int DOF = 6;

int main()
{
	srand(time(NULL));

	NdofOpenChain< DOF > robot;
	StatePtr state = robot.makeState();
	//renderer.getViewer().realize();

	//for (unsigned int i = 0; i < DOF; i++)
	//{
	//	state->setJointStatePos(i, (rand()%100 - 50)/100.0 * PI);
	//	state->setJointStateVel(i, (rand() % 100 - 50) / 100.0 * PI);
	//}
	{
		//state->setJointStatePos(0, 1);
		//state->setJointStatePos(1, 2);
		//state->setJointStatePos(2, PI_HALF);

		for (unsigned int i = 0; i < DOF; i++)
		{
			state->setJointStatePos(i, (rand()%100 - 50)/100.0 * PI_HALF);
		}

		state->setJointStateVel(0, 1);
		state->setJointStateVel(1, 2);
		state->setJointStateVel(2, 3);

		state->setJointStateAcc(0, 1);
		state->setJointStateAcc(1, 2);
		state->setJointStateAcc(2, 3);

		robot.solveForwardKinematics(*state);
		robot.solveDiffForwardKinematics(*state);
		robot.solve2ndDiffForwardKinematics(*state);

		robot.solveJacobian(*state);
		robot.solveJacobianDot(*state);

		for (unsigned int i = 0; i < DOF; i++)
		{
			cout << "q(" << i << ") = " << state->getJointStatePos(i) << endl;
		}
		for (unsigned int i = 0; i < DOF; i++)
		{
			cout << "dq(" << i << ") = " << state->getJointStateVel(i) << endl;
		}
		for (unsigned int i = 0; i < DOF; i++)
		{
			cout << "ddq(" << i << ") = " << state->getJointStateAcc(i) << endl;
		}

		cout << "TRANSFORM: " << endl << state->getLinkStateSE3(3) << endl << endl;
		cout << "VELOCITY: " << endl << state->getLinkStateVel(3) << endl << endl;
		cout << "ACCELERATION: " << endl << state->getLinkStateAcc(3) << endl << endl;

		cout << "JACOBIAN: " << endl;
		for (unsigned int i = 0; i < DOF; i++)
		{
			cout << state->getJointStateScrew(i) << endl;
		}
		cout << endl;

		cout << "JACOBIAN_DOT: " << endl;
		for (unsigned int i = 0; i < DOF; i++)
		{
			cout << state->getJointStateScrewDot(i) << endl;
		}
		cout << endl;
	}
	PERFORM_TEST(
		for (unsigned int i = 0; i < DOF; i++)
		{
			state->setJointStatePos(i, (rand() % 100 - 50) / 100.0 * PI);
		};
	robot.solveForwardKinematics(*state);,
		10e+6
		);

	//OSG_simpleRender renderer(robot, *state, 600, 600);
	//renderer.getViewer().run();

	_getch();

	return 0;
}