#include <iostream>
#include <fstream>
#include <conio.h>
#include <rovin\TimeOptimization\TOPP.h>

#include "efortRobot.h"

using namespace std;
using namespace rovin;

int main()
{	
	SerialOpenChainPtr robot(new efortRobot());
	StatePtr state = robot->makeState();

	// data 
	ifstream input("trajectory.txt");
	if (input.fail()) cout << "파일 열기 실패" << endl;
	else cout << "파일 열기 성공" << endl;

	Real tmp;
	unsigned int cnt = 0;
	while (!input.eof()) {
		input >> tmp;
		cnt++;
	}
	input.close();

	ifstream trajectory("trajectory.txt");
	if (trajectory.fail()) cout << "파일 열기 실패" << endl;
	else cout << "파일 열기 성공" << endl;

	unsigned int dof = 6;
	unsigned int data_num = cnt / dof;
	MatrixX q(dof, data_num);

	for (int i = 0; i < data_num; i++)
		for (int j = 0; j < dof; j++)
			trajectory >> q(j,i);

	cout << q.col(0) << endl;
	cout << q.col(q.cols()-1) << endl;
	
	trajectory.close();

	//MatrixX q_data(dof, 7);
	//for (int i = 0; i < 7; i++)
	//	for (int j = 0; j < dof; j++)
	//		q_data(j, i) = q(j, 80*i);

	// Time optimization
	Real ds = 1e-3;
	Real vi = 0;
	Real vf = 0;
	Real si = 0;
	Real sf = 1;
		
	//TOPP topp(q_data, robot, ds, vi, vf, si, sf);
	TOPP topp(q, robot, ds, vi, vf, si, sf);

	//topp.saveMVCandSP2txt();

//	topp.saveMVCandSP2txt();

	// generate trajectory
	topp.generateTrajectory();

	cout << "Final time : " << topp.getFinalTime() << endl;



	_getch();
	return 0;
}