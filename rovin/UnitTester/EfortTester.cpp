#include <cstdio>
#include <conio.h>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include <rovin/Renderer/OSG_simpleRender.h>

#include "efortRobot.h"
#include <rovin/Utils/Diagnostic.h>


using namespace std;
using namespace rovin;

int main()
{
	//efortRobot robot;
	//rovin::StatePtr state = robot.makeState();

	//const int ParameterN = 3;


	//MatrixX dqdp(6, ParameterN), dqdotdp(6, ParameterN), dqddotdp(6, ParameterN);
	//dqdp.setRandom(); dqdotdp.setRandom(); dqddotdp.setRandom();

	/*FILE *out;
	fopen_s(&out, "C:\\Users\\근준\\Documents\\rovin\\rovin\\data.txt", "wt");
	for (int i = 0; i < 6; i++) fprintf_s(out, "%lf\n", q(i));
	for (int i = 0; i < 6; i++) fprintf_s(out, "%lf\n", qdot(i));
	for (int i = 0; i < 6; i++) fprintf_s(out, "%lf\n", qddot(i));
	for (int i = 0; i < 6; i++) for (int j = 0; j < ParameterN; j++) fprintf_s(out, "%lf\n", dqdp(i, j));
	for (int i = 0; i < 6; i++) for (int j = 0; j < ParameterN; j++) fprintf_s(out, "%lf\n", dqdotdp(i, j));
	for (int i = 0; i < 6; i++) for (int j = 0; j < ParameterN; j++) fprintf_s(out, "%lf\n", dqddotdp(i, j));
	fcloseall();*/

	//robot.solveInverseDynamics(*state);
	//cout << state->getJointStateTorque() << endl;
	//cout << robot.solveDiffInverseDynamics(*state, dqdp, dqdotdp, dqddotdp) << endl;

	//robot.solveForwardKinematics(*state);
	//OSG_simpleRender renderer(robot, *state, 600, 600);
	//renderer.getViewer().run();

	//VectorX q(6), qdot(6), qddot(6);
	//state->setJointStatePos(q);
	//state->setJointStateVel(qdot);
	//state->setJointStateAcc(qddot);


	// file input output
	ifstream input("trajectory.txt");
	if (input.fail()) cout << "파일 열기 실패" << endl;
	else cout << "파일 열기 성공" << endl;

	float tmp;
	unsigned int cnt = 0;
	while (!input.eof()) {
		input >> tmp;
		cnt++;
	}
	cout << cnt << endl;
	input.close();

	ifstream trajectory("trajectory.txt");
	if (trajectory.fail()) cout << "파일 열기 실패" << endl;
	else cout << "파일 열기 성공" << endl;

	unsigned int dof = 6;
	unsigned int data_num = cnt / dof;

	std::vector<VectorX> q(dof);
	for (int i = 0; i < dof; i++)
		q[i].resize(data_num);

	for (int i = 0; i < data_num; i++)
		for (int j = 0; j < dof; j++)
			trajectory >> q[j](i);

	trajectory.close();


	// rendering
	efortRobot robot;
	rovin::StatePtr state = robot.makeState();

	state->setJointStatePos(0, q[0](0));
	state->setJointStatePos(1, q[1](0));
	state->setJointStatePos(2, q[2](0));
	state->setJointStatePos(3, q[3](0));
	state->setJointStatePos(4, q[4](0));
	state->setJointStatePos(5, q[5](0));
	robot.solveForwardKinematics(*state);

	OSG_simpleRender renderer(robot, *state, 600, 600);
	renderer.getViewer().realize();

	double frameRate = 50;

	cnt = 0;
	double c = clock();
	while (1)
	{
		if (clock() - c >= 1000 / frameRate)
		{
			if (cnt == data_num) cnt = 0;
			state->setJointStatePos(0, q[0](cnt));
			state->setJointStatePos(1, q[1](cnt));
			state->setJointStatePos(2, q[2](cnt));
			state->setJointStatePos(3, q[3](cnt));
			state->setJointStatePos(4, q[4](cnt));
			state->setJointStatePos(5, q[5](cnt));
			cnt++;
			robot.solveForwardKinematics(*state);
		}
		renderer.updateFrame();
	}


	_getch();

	return 0;
}