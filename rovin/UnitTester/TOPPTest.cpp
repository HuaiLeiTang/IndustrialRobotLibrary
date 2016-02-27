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
	loadData(q_data);

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

//void GivenPathTimeOptimization::InputGivenPath(const std::vector<SE3, Eigen::aligned_allocator<SE3>>& givenPath)
//{
//	_givenPath = givenPath;
//
//	int pathsize = givenPath.size();
//	_s.resize(pathsize);
//	for (int i = 0; i < pathsize; i++)
//		_s[i] = Real(i) / Real(pathsize - 1);
//
//	_sdot.resize(pathsize);
//	_sddot.resize(pathsize);
//
//	solveInvKinAll();

	//int pathsize = _givenPath.size();
	//LOGIF((pathsize != 0), "GivenPathTimeOptimization::solveInvKinAll error : No given path trajectory");

	//StatePtr statePtr = _socRobotPtr->makeState();
	//int dof = statePtr->getDof();

	/////< calculate qs: inverse kinematics at each given location(SE3)
	//_q.resize(pathsize);
	//for (int i = 0; i < pathsize; i++)
	//{
	//	_socRobotPtr->solveInverseKinematics(*statePtr, _givenPath[i]);
	//	_q[i] = statePtr->getJointStatePos();
	//}
//}
