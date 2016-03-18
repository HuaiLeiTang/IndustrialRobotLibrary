#include <iostream>
#include <fstream>
#include <conio.h>
#include <rovin\TimeOptimization\TOPP.h>
#include <rovin\TimeOptimization\P2PTimeOptTemp.h>
#include <rovin\Renderer\OSG_simpleRender.h>

#include "efortRobot.h"
#include <string>

using namespace std;
using namespace rovin;

void loadData(MatrixX& data);

SerialOpenChainPtr robot(new efortRobot());
StatePtr state;
unsigned int dof;

int main()
{	
	state = robot->makeState();
	dof = robot->getNumOfJoint();

	//////////////////////////////////////////////////////////////////////
	// AVP test
	AVP_RRT avp_rrt(robot, CONSTRAINT_TYPE::TORQUE);

	// start and goal
	VectorX qs(dof), qg(dof), qsdot(dof), qgdot(dof);
	qs << 0.115499, 0.220374, 0.151515, 0.15633, 0.0438677, 0.0182414;
	qg << 1.15739, 1.11361, 1.04554, 1.25135, 1.0009, 0.790497;
	qsdot.setZero();
	qgdot.setZero();

	VectorX q1(dof), q1dot(dof);
	q1 << 0.638163, 0.604759, 0.800566, 0.846903, 0.488474, 0.422513;
	q1dot.setZero();

	std::vector<WayPoint> wayPt(3);
	wayPt[0].setJointq(qs);
	wayPt[0].setJointqdot(qsdot);
	wayPt[1].setJointq(q1);
	wayPt[1].setJointqdot(q1dot);
	wayPt[wayPt.size() - 1].setJointq(qg);
	wayPt[wayPt.size() - 1].setJointqdot(qgdot);

	avp_rrt.setWayPoints(wayPt);
	avp_rrt.generateTrajectory();

	MatrixX finalTraj = avp_rrt.getFinalPath();

	//rendering
	//int num1 = 100;
	//finalTraj.resize(6, 3* num1);
	//for (int i = 0; i < num1; i++)
	//{
	//	finalTraj.col(num1*0 + i) = qs;
	//	finalTraj.col(num1*1 + i) = q1;
	//	finalTraj.col(num1*2 + i) = qg;
	//}

	state->setJointStatePos(0, finalTraj(0, 0));
	state->setJointStatePos(1, finalTraj(1, 0));
	state->setJointStatePos(2, finalTraj(2, 0));
	state->setJointStatePos(3, finalTraj(3, 0));
	state->setJointStatePos(4, finalTraj(4, 0));
	state->setJointStatePos(5, finalTraj(5, 0));
	robot->solveForwardKinematics(*state);

	OSG_simpleRender renderer(*robot, *state, 600, 600);
	renderer.getViewer().realize();

	double frameRate = 50;

	int cnt = 0;
	int data_num = finalTraj.cols();
	double c = clock();
	while (1)
	{
		if (clock() - c >= 1000 / frameRate)
		{
			if (cnt == data_num) cnt = 0;
			state->setJointStatePos(0, finalTraj(0, cnt));
			state->setJointStatePos(1, finalTraj(1, cnt));
			state->setJointStatePos(2, finalTraj(2, cnt));
			state->setJointStatePos(3, finalTraj(3, cnt));
			state->setJointStatePos(4, finalTraj(4, cnt));
			state->setJointStatePos(5, finalTraj(5, cnt));
			cnt++;
			robot->solveForwardKinematics(*state);
		}
		renderer.updateFrame();
	}

	cout << "Program complete" << endl;
	_getch();
	return 0;
}

void loadData(MatrixX& data)
{
	ifstream input("D:/jkkim/Documents/trajectory_wy.txt");
	if (input.fail()) cout << "파일 열기 실패" << endl;
	else cout << "파일 열기 성공" << endl;

	Real tmp;
	unsigned int cnt = 0;
	while (!input.eof()) {
		input >> tmp;
		cnt++;
	}
	input.close();

	ifstream trajectory("D:/jkkim/Documents/trajectory_wy.txt");
	if (trajectory.fail()) cout << "파일 열기 실패" << endl;
	else cout << "파일 열기 성공" << endl;

	unsigned int dof = 6;
	unsigned int data_num = cnt / dof;
	MatrixX q(dof, data_num);

	for (unsigned int i = 0; i < data_num; i++)
		for (unsigned int j = 0; j < dof; j++)
			trajectory >> q(j, i);

	trajectory.close();

	data = q;
}