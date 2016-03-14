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
	////////////////////////////////////////////////////////////////////////
	// TOPP test

	//std::list<Real> a;
	//std::list<Real> b;

	//a.push_back(10);
	//a.push_back(20);
	//a.push_front(3);
	//b.push_front(-3);
	//b.push_back(-19);
	//b.push_front(-10);

 //   a.insert(a.begin(), b.begin(), b.end());
	//b.reverse();

	//std::list<Real>::iterator a_iter = a.begin();
	//std::list<Real>::iterator b_iter = b.begin();

	//for (a_iter = a.begin(); a_iter != a.end(); a_iter++)
	//	std::cout << (*a_iter) << '\t';
	//std::cout << std::endl;
	//for (b_iter = b.begin(); b_iter != b.end(); b_iter++)
	//	std::cout << (*b_iter) << '\t';

	//state = robot->makeState();
	//dof = robot->getNumOfJoint();

	//MatrixX q_data;
	//loadData(q_data);

	//// Time optimization
	//Real ds = 1e-3, vi = 0, vf = 0, si = 0, sf = 1;

	//// consider only torque constraint
	//TOPP topp1(q_data, robot, vi, vf, ds, si, sf, CONSTRAINT_TYPE::TORQUE_VEL_ACC);
	//topp1.generateTrajectory();
	//cout << "[ consider only torque constraint ]" << endl;
	//cout << "Final time : " << topp1.getFinalTime() << endl << endl;

	//// consider torque, velocity and acceleration constraint
	//TOPP topp2(q_data, robot, vi, vf, ds, si, sf, CONSTRAINT_TYPE::TORQUE_VEL_ACC);
	//topp2.generateTrajectory();
	//cout << "[ consider torque, velocity and acceleration constraint ]" << endl;
	//cout << "Final time : " << topp2.getFinalTime() << endl << endl;

	////////////////////////////////////////////////////////////////////////
	// AVP test
	//AVP_RRT avp_rrt(robot, CONSTRAINT_TYPE::TORQUE_VEL_ACC);

	// start and goal
	//VectorX qs(dof), qg(dof), qsdot(dof), qgdot(dof);
	//qs << 0.115499, 0.220374, 0.151515, 0.15633, 0.0438677, 0.0182414;
	//qg << 1.15739, 1.11361, 1.04554, 1.25135, 1.0009, 0.790497;
	//qsdot.setZero();
	//qgdot.setZero();

	//VectorX q1(dof), q1dot(dof);
	//q1 << 0.638163, 0.604759, 0.800566, 0.846903, 0.488474, 0.422513;
	//q1dot.setZero();

	//std::vector<WayPoint> wayPt(3);
	//wayPt[0].setJointq(qs);
	//wayPt[0].setJointqdot(qsdot);
	//wayPt[1].setJointq(q1);
	//wayPt[1].setJointqdot(q1dot);
	//wayPt[wayPt.size() - 1].setJointq(qg);
	//wayPt[wayPt.size() - 1].setJointqdot(qgdot);

	//AVP_RRT avp_rrt(robot, CONSTRAINT_TYPE::TORQUE_VEL_ACC);
	//avp_rrt.setWayPoints(wayPt);
	//avp_rrt.generateTrajectory();

	//MatrixX finalTraj = avp_rrt.getFinalPath();
	//int data_num = finalTraj.cols();

	//MatrixX q_data;
	//loadData(q_data);
	//std::list<VectorX> Pnew;
	//for (int i = 0; i < q_data.cols(); i++)
	//	Pnew.push_back(q_data.col(i));

	//Vector2 nearInterval(0, 0.5);
	//Vector2 endInterval;
	//avp_rrt.runAVP(Pnew, nearInterval, endInterval);
	////avp_rrt.runAVPbackward(Pnew, nearInterval, endInterval);

	////////////////////////////////////////////////////////////////////////
	// Interpolation test
	MatrixX q_data;
	loadData(q_data);
	VectorX qnear = q_data.col(0);
	VectorX qnearVel = VectorX::Ones(6) * 0.1;
	VectorX qrand = q_data.col(30);
	Real dist = 1.0;

	cout << "[qnear]" << endl;
	cout << qnear << endl << endl;
	cout << "[qrand]" << endl;
	cout << qrand << endl << endl;
	cout << "[qnearVel]" << endl;
	cout << qnearVel << endl << endl;


	std::list<VectorX> Pnew;
	VectorX qnew;

	Vertex* nVertex = new Vertex();
	nVertex->setconfig(qnear);
	nVertex->setconfigVel(qnearVel);

	//Vector2 nearInterval(0, 1);
	//Vector2 endInterval;
	//avp_rrt.runAVP(Pnew, nearInterval, endInterval);


	// rendering
	//std::cout << "---";
	//int num1 = 100;
	//finalTraj.resize(6, 3* num1);
	//for (int i = 0; i < num1; i++)
	//{
	//	finalTraj.col(num1*0 + i) = qs;
	//	finalTraj.col(num1*1 + i) = q1;
	//	finalTraj.col(num1*2 + i) = qg;
	//}

	//state->setJointStatePos(0, finalTraj(0, 0));
	//state->setJointStatePos(1, finalTraj(1, 0));
	//state->setJointStatePos(2, finalTraj(2, 0));
	//state->setJointStatePos(3, finalTraj(3, 0));
	//state->setJointStatePos(4, finalTraj(4, 0));
	//state->setJointStatePos(5, finalTraj(5, 0));
	//robot->solveForwardKinematics(*state);

	//OSG_simpleRender renderer(*robot, *state, 600, 600);
	//renderer.getViewer().realize();

	//double frameRate = 50;

	//int cnt = 0;
	//double c = clock();
	//while (1)
	//{
	//	if (clock() - c >= 1000 / frameRate)
	//	{
	//		if (cnt == data_num) cnt = 0;
	//		state->setJointStatePos(0, finalTraj(0, cnt));
	//		state->setJointStatePos(1, finalTraj(1, cnt));
	//		state->setJointStatePos(2, finalTraj(2, cnt));
	//		state->setJointStatePos(3, finalTraj(3, cnt));
	//		state->setJointStatePos(4, finalTraj(4, cnt));
	//		state->setJointStatePos(5, finalTraj(5, cnt));
	//		cnt++;
	//		robot->solveForwardKinematics(*state);
	//	}
	//	renderer.updateFrame();
	//}

	AVP_RRT avp_rrt(robot, CONSTRAINT_TYPE::TORQUE_VEL_ACC);
	avp_rrt.interpolate(nVertex, qrand, dist, Pnew, qnew);

	cout << "[qnew]" << endl;
	cout << qnew << endl << endl;

	delete nVertex;
	
	cout << "Program complete" << endl;
	_getch();
	return 0;
}

void loadData(MatrixX& data)
{
	//ifstream input("C:/Users/ksh/Documents/trajectory_wy.txt");
	//ifstream input("D:/jkkim/Documents/trajectory_wy.txt");
	ifstream input("C:/Users/crazy/Desktop/Time optimization/trajectory text/trajectory_wy.txt");
	if (input.fail()) cout << "파일 열기 실패" << endl;
	else cout << "파일 열기 성공" << endl;

	Real tmp;
	unsigned int cnt = 0;
	while (!input.eof()) {
		input >> tmp;
		cnt++;
	}
	input.close();

	//ifstream trajectory("C:/Users/ksh/Documents/trajectory_wy.txt");
	//ifstream trajectory("D:/jkkim/Documents/trajectory_wy.txt");
	ifstream trajectory("C:/Users/crazy/Desktop/Time optimization/trajectory text/trajectory_wy.txt");
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