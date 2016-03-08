#include <iostream>
#include <fstream>
#include <conio.h>
#include <rovin\TimeOptimization\TOPP.h>
#include <rovin\TimeOptimization\P2PTimeOptimization.h>

#include "efortRobot.h"
#include <string>

using namespace std;
using namespace rovin;

int main()
{
	SerialOpenChainPtr robot(new efortRobot());
	StatePtr state;
	unsigned int dof;

	state = robot->makeState();
	dof = robot->getNumOfJoint();









	_getch();
	return 0;
}
