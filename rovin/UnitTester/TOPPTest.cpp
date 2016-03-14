#include <iostream>
#include <fstream>
#include <conio.h>
#include <rovin\TimeOptimization\TOPP.h>
#include <rovin\TimeOptimization\P2PTimeOptTemp.h>

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