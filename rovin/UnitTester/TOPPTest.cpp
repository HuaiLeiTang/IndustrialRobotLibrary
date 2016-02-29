#include <iostream>
#include <fstream>
#include <conio.h>
#include <rovin\TimeOptimization\TOPP.h>

#include "efortRobot.h"

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

	// generate trajectory
	MatrixX q_data;
	loadData(q_data);

	// rendering


	// Time optimization
	Real ds = 1e-3, vi = 0, vf = 0, si = 0, sf = 1;
	TOPP topp(q_data, robot, ds, vi, vf, si, sf);
	//topp.generateTrajectory();
	topp.saveMVCandSP2txt();

	_getch();
	return 0;
}

void loadData(MatrixX& data)
{
	//ifstream input("trajectory.txt");
	//ifstream input("C:/Users/crazy/Desktop/Time optimization/trajectory text/trajectory_wy.txt");
	ifstream input("C:/Users/ksh/Documents/trajectory_wy.txt");
	if (input.fail()) cout << "파일 열기 실패" << endl;
	else cout << "파일 열기 성공" << endl;

	Real tmp;
	unsigned int cnt = 0;
	while (!input.eof()) {
		input >> tmp;
		cnt++;
	}
	input.close();

	//ifstream trajectory("trajectory.txt");
	//ifstream trajectory("C:/Users/crazy/Desktop/Time optimization/trajectory text/trajectory_wy.txt");
	ifstream trajectory("C:/Users/ksh/Documents/trajectory_wy.txt");
	if (trajectory.fail()) cout << "파일 열기 실패" << endl;
	else cout << "파일 열기 성공" << endl;

	unsigned int dof = 6;
	unsigned int data_num = cnt / dof;
	MatrixX q(dof, data_num);

	for (int i = 0; i < data_num; i++)
		for (int j = 0; j < dof; j++)
			trajectory >> q(j, i);

	trajectory.close();

	data = q;

	//data = MatrixX(dof, 7);
	//for (int i = 0; i < 7; i++)
	//	for (int j = 0; j < dof; j++)
	//		data(j, i) = q(j, 80 * i);
}
