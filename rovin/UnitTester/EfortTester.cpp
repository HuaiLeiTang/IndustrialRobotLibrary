#include <cstdio>
#include <conio.h>
#include <ctime>`
#include <cstdlib>
#include <iostream>
#include <rovin/Renderer/OSG_simpleRender.h>

#include "efortRobot.h"
#include <rovin/Utils/Diagnostic.h>

#include <rovin/TimeOptimization/GivenPathTimeOptimization.h>

using namespace std;
using namespace rovin;

int main()
{
	SerialOpenChainPtr robot = SerialOpenChainPtr(new efortRobot);
	rovin::StatePtr state = robot->makeState();

	vector<SE3, Eigen::aligned_allocator<SE3>> T_traj;
	SE3 T_temp;

	FILE *in;
	fopen_s(&in, "SE3.txt", "r");
	while (!feof(in))
	{
		MatrixX T(3, 4);
		fscanf_s(in, "%lf", &T(0, 0));
		if (feof(in)) break;

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				if (i == 0 && j == 0) continue;
				fscanf_s(in, "%lf", &T(i, j));
			}
		}
		T_temp.setRotation(SO3::Projection(T.block<3, 3>(0, 0)));
		T_temp.setPosition(T.block<3, 1>(0, 3));
		T_traj.push_back(T_temp);
	}

	GivenPathTimeOptimization GPTO(robot, T_traj);

	GPTO.solveMinimumTimeOptimization(0.0, 0.0);

	_getch();

	//OSG_simpleRender renderer(*robot, *state, 600, 600);
	//renderer.getViewer().realize();

	//int count = 0;
	//int t = clock();
	//while(1)
	//{
	//	if (clock() - t >= 30)
	//	{
	//		t = clock();
	//		state->setJointStatePos(GPTO.getq(count));
	//		robot->solveForwardKinematics(*state);
	//		count++;
	//		if (count >= T_traj.size()) count = 0;
	//	}
	//	renderer.updateFrame();
	//}

	return 0;
}