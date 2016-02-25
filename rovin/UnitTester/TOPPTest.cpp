#include <iostream>
#include <conio.h>
#include <rovin\TimeOptimization\TOPP.h>

#include "efortRobot.h"

using namespace std;
using namespace rovin;

int main()
{	
	efortRobot robot;
	StatePtr state = robot.makeState();



	_getch();
	return 0;
}