#include <cstdio>
#include <conio.h>
#include <ctime>`
#include <cstdlib>
#include <iostream>
#include <rovin/Renderer/OSG_simpleRender.h>

#include "efortRobot.h"
#include <rovin/Utils/Diagnostic.h>

#include <cstdio>

using namespace std;
using namespace rovin;

int main()
{
	// file input
	vector<VectorX> q_result_vec;
	VectorX q_temp(6);

	FILE *in;
	fopen_s(&in, "q_result.txt", "r");
	while (!feof(in))
	{
		for (int i = 0; i < 6; i++)
		{
			fscanf_s(in, "%lf", &q_temp(i));
		}
		q_result_vec.push_back(q_temp);
	}


	// robot animation
	SerialOpenChainPtr robot = SerialOpenChainPtr(new efortRobot);
	rovin::StatePtr state = robot->makeState();

	OSG_simpleRender renderer(*robot, *state, 600, 600);
	renderer.getViewer().realize();

	int count = 0;
	int t = clock();
	while(1)
	{
		if (clock() - t >= 10)
		{
			t = clock();
			state->setJointStatePos(q_result_vec[0]);
			robot->solveForwardKinematics(*state);
			count++;
			if (count >= q_result_vec.size())
				count = 0;
		}
		renderer.updateFrame();
	}

	//const int ParameterN = 3;

	//VectorX q(6), qdot(6), qddot(6);
	//MatrixX dqdp(6, ParameterN), dqdotdp(6, ParameterN), dqddotdp(6, ParameterN);

	//q.setRandom(); qdot.setRandom(); qddot.setRandom();
	//dqdp.setRandom(); dqdotdp.setRandom(); dqddotdp.setRandom();

	//FILE *out;
	//fopen_s(&out, "C:\\Users\\±Ÿ¡ÿ\\Documents\\rovin\\rovin\\data.txt", "wt");
	//for (int i = 0; i < 6; i++) fprintf_s(out, "%lf\n", q(i));
	//for (int i = 0; i < 6; i++) fprintf_s(out, "%lf\n", qdot(i));
	//for (int i = 0; i < 6; i++) fprintf_s(out, "%lf\n", qddot(i));
	//for (int i = 0; i < 6; i++) for (int j = 0; j < ParameterN; j++) fprintf_s(out, "%lf\n", dqdp(i, j));
	//for (int i = 0; i < 6; i++) for (int j = 0; j < ParameterN; j++) fprintf_s(out, "%lf\n", dqdotdp(i, j));
	//for (int i = 0; i < 6; i++) for (int j = 0; j < ParameterN; j++) fprintf_s(out, "%lf\n", dqddotdp(i, j));
	//fcloseall();

	//state->setJointStatePos(q);
	//state->setJointStateVel(qdot);
	//state->setJointStateAcc(qddot);

	//robot.solveInverseDynamics(*state);
	//cout << state->getJointStateTorque() << endl;
	//cout << robot.differentiateInverseDynamics(*state, dqdp, dqdotdp, dqddotdp) << endl;

	//robot.solveForwardKinematics(*state);

	//OSG_simpleRender renderer(robot, *state, 600, 600);
	//renderer.getViewer().run();

	_getch();

	return 0;
}