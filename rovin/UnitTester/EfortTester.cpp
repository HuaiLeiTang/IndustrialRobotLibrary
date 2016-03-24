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
	/*
	* FILE INPUT
	*/
	// Matrices to save file input
	vector<VectorX> toolpath_vec;
	vector<VectorX> q_result_vec;
	vector<VectorX> q_ref_result_vec;
	VectorX VecX_temp(6);


	vector<Real> time_result_vec;
	Real time_temp;

	FILE *in;
	// tool path
	fopen_s(&in, "selectedPath.txt", "r");
	while (!feof(in))
	{
		for (int i = 0; i < 6; i++)
		{
			fscanf_s(in, "%lf", &VecX_temp(i));
		}
		toolpath_vec.push_back(VecX_temp * 0.001); // scaling
	}

	// optimal joint trajectory
	fopen_s(&in, "q_result.txt", "r");
	while (!feof(in))
	{
		for (int i = 0; i < 6; i++)
		{
			fscanf_s(in, "%lf", &VecX_temp(i));
		}
		q_result_vec.push_back(VecX_temp);
	}

	// reference joint trajectory
	fopen_s(&in, "q_ref_result.txt", "r");
	while (!feof(in))
	{
		for (int i = 0; i < 6; i++)
		{
			fscanf_s(in, "%lf", &VecX_temp(i));
		}
		q_ref_result_vec.push_back(VecX_temp);
	}

	// time array
	fopen_s(&in, "time_result.txt", "r");
	while (!feof(in))
	{
		fscanf_s(in, "%lf", &time_temp);
		time_result_vec.push_back(time_temp);
	}
	fcloseall();



	/*
	* RENDERING
	*/
	// Renderer setting
	SerialOpenChainPtr robot = SerialOpenChainPtr(new efortRobot);
	rovin::StatePtr state = robot->makeState();
	rovin::StatePtr state2 = robot->makeState();

	OSG_simpleRender renderer(*robot, *state, 600, 600);
	renderer.getViewer().realize();

	OSG_simpleRender renderer2(*robot, *state2, 600, 600);
	renderer2.getViewer().realize();

	// Render trajectory
	Line	EE_pathLine;
	renderer.addGeometry(EE_pathLine);
	renderer2.addGeometry(EE_pathLine);
	//EE_pathLine.setColor(0 / 255.0, 0 / 255.0, 255 / 255.0, 1.0);

	Vector3 defaultPos;
	defaultPos << -0.500, 0, 0;
	for (int i = 0; i < toolpath_vec.size(); i++)
	{
		osg::Vec3 pos_temp;
		for (int j = 0; j < 3; j++)
			pos_temp[j] = (toolpath_vec[i])[j] + defaultPos[j];
		EE_pathLine.push_back(pos_temp);
		//for (int ii = 0; ii < 3; ii++)
		//	cout << pos_temp[ii] << "   ";
		//cout << endl;
	}


	// Render robot animation
	int count = 0;
	int t0 = clock();
	Real speed_render = 1.0;
	while(1)
	{
		if (speed_render * Real(clock() - t0) / CLOCKS_PER_SEC >= time_result_vec[count])
		{
			state->setJointStatePos(q_result_vec[count]);
			robot->solveForwardKinematics(*state);

			state2->setJointStatePos(q_ref_result_vec[count]);
			robot->solveForwardKinematics(*state2);

			count++;
			if (count >= q_result_vec.size())
			{
				count = 0;
				t0 = clock();
			}
		}
		renderer.updateFrame();
		renderer2.updateFrame();
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