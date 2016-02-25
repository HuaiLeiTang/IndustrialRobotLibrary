#include <iostream>
#include <fstream>
#include <conio.h>
#include <rovin\TimeOptimization\TOPP.h>

#include "efortRobot.h"

using namespace std;
using namespace rovin;

int main()
{	
	efortRobot robot;
	StatePtr state = robot.makeState();

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
	MatrixX q(data_num, dof);

	for (int i = 0; i < data_num; i++)
		for (int j = 0; j < dof; j++)
			trajectory >> q(i,j);
	
	trajectory.close();

	MatrixX q_data(6, dof);
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < dof; j++)
		{
			q_data(i, j) = q(50*i, j);
		}
	}


	_getch();
	return 0;
}