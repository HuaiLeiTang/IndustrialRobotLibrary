#include <cstdio>
#include <conio.h>
#include <ctime>`
#include <cstdlib>
#include <iostream>
#include <rovin/Renderer/OSG_simpleRender.h>

#include "efortRobot.h"
#include <rovin/Utils/Diagnostic.h>

using namespace std;
using namespace rovin;

int main()
{
	efortRobot robot;
	rovin::StatePtr state = robot.makeState();

	VectorX q(6);
	q.setZero();
	state->setJointStatePos(q);

	robot.solveForwardKinematics(*state);

	//OSG_simpleRender renderer(robot, *state, 600, 600);
	//renderer.getViewer().run();

	return 0;
}