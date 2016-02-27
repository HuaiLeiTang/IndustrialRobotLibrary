#include <iostream>
#include <fstream>
#include <conio.h>
#include <rovin\TimeOptimization\TOPP.h>

#include "efortRobot.h"

using namespace std;
using namespace rovin;

void loadData(MatrixX& data);
void generateTrajectory(MatrixX& data);

SerialOpenChainPtr robot(new efortRobot());
StatePtr state;
unsigned int dof;

int main()
{	
	state = robot->makeState();
	dof = robot->getNumOfJoint();

	// generate trajectory
	MatrixX q_data;


	// rendering

	
	// Time optimization
	Real ds = 1e-3, vi = 0, vf = 0, si = 0, sf = 1;
	TOPP topp(q_data, robot, ds, vi, vf, si, sf);
	topp.generateTrajectory();
	cout << "Final time : " << topp.getFinalTime() << endl;

	_getch();
	return 0;
}





void loadData(MatrixX& data)
{
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
			trajectory >> q(j, i);

	trajectory.close();

	data = MatrixX(dof, 7);
	for (int i = 0; i < 7; i++)
		for (int j = 0; j < dof; j++)
			data(j, i) = q(j, 80 * i);
}

void generateTrajectory(MatrixX& data)
{
	VectorX q(dof), qdot(dof), qddot(dof);
	SE3 GaolT;

	q << 1.5, 0.8, 1.0, -1.5, 0.1, 1.4;
	qdot.setZero(); qddot.setZero();
	state->setJointStatePos(q);
	state->setJointStateVel(qdot);
	state->setJointStateAcc(qddot);



}