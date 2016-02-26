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
	MatrixX q(data_num, dof);

	for (int i = 0; i < data_num; i++)
		for (int j = 0; j < dof; j++)
			trajectory >> q(i,j);
	
	trajectory.close();

	MatrixX q_data(dof, 7);
	for (int i = 0; i < 7; i++)
		for (int j = 0; j < dof; j++)
			q_data(j, i) = q(80*i, j);

	// Time optimization
	Real ds = 1e-3;
	Real vi = 0;
	Real vf = 0;
	Real si = 0;
	Real sf = 1;

	TOPP topp(q_data, robot, ds, vi, vf, si, sf);

	// cout << topp._torqueConstraint << endl;

	//topp.saveMVCandSP2txt();

	// generate trajectory
	topp.generateTrajectory();
	

	/////////////////////////////////////////////////////// experiment
	//cout << "determintAlphaBeta" << endl;
	//Vector2 result1 = topp.determineAlphaBeta(0.23, 2.795);

	//cout << endl;
	//cout << "calculateBandC" << endl;
	//std::vector<VectorX> result2 = topp.calculateBandC(0.23);

	//cout << "result2[0]" << endl;
	//cout << "result2[0] : " << result2[0].size() << endl;
	//cout << result2[0] << endl;
	//cout << "result2[1]" << endl;
	//cout << "result2[1] : " << result2[1].size() << endl;
	//cout << result2[1] << endl;

	//cout << "result2[0] * 2.795*2.795 + result2[1]" << endl;
	//cout << -result2[0] * 2.795*2.795 - result2[1] << endl;



	//cout << endl;
	//cout << "calculateMVCPoint" << endl;
	//Real result2 = topp.calculateMVCPoint(0.23);
	
	//cout << topp._tf_result << endl;





	//////////////////////////////////////////
	//VectorX q(6), dq(6), ddq(6);
	//q.setOnes();
	//dq.setOnes();
	//ddq.setOnes();
	//state->setJointStatePos(q);
	//state->setJointStateVel(dq);
	//state->setJointStateAcc(ddq);
	//robot->solveInverseDynamics(*state);
	//VectorX MCg = state->getJointStateTorque();
	//
	//dq.setZero(); ddq.setZero();
	//state->setJointStatePos(q);
	//state->setJointStateVel(dq);
	//state->setJointStateAcc(ddq);
	//robot->solveInverseDynamics(*state);
	//VectorX g = state->getJointStateTorque();

	//dq.setOnes();
	//state->setJointStatePos(q);
	//state->setJointStateVel(dq);
	//state->setJointStateAcc(ddq);
	//robot->solveInverseDynamics(*state);
	//VectorX Cg = state->getJointStateTorque();

	//VectorX C = Cg - g;
	//VectorX MC = MCg - g;

	//Real s = 1.5;
	//q.setOnes();
	//dq.setOnes(); dq = dq * s;
	//ddq.setOnes(); ddq = ddq * s *s;
	//state->setJointStatePos(q);
	//state->setJointStateVel(dq);
	//state->setJointStateAcc(ddq);
	//robot->solveInverseDynamics(*state);
	//VectorX _MCg = state->getJointStateTorque();

	//dq.setZero(); ddq.setZero();
	//state->setJointStatePos(q);
	//state->setJointStateVel(dq);
	//state->setJointStateAcc(ddq);
	//robot->solveInverseDynamics(*state);
	//VectorX _g = state->getJointStateTorque();

	//dq.setOnes();
	//dq = dq*s;
	//state->setJointStatePos(q);
	//state->setJointStateVel(dq);
	//state->setJointStateAcc(ddq);
	//robot->solveInverseDynamics(*state);
	//VectorX _Cg = state->getJointStateTorque();

	//VectorX _C = _Cg - _g;
	//VectorX _MC = _MCg - _g;
	//
	//cout << "MC" << endl;
	//cout << MC*s*s << endl;
	//cout << "_MC" << endl;
	//cout << _MC << endl;

	//cout << "g" << endl;
	//cout << g << endl;
	//cout << "_g" << endl;
	//cout << _g << endl;

	_getch();
	return 0;
}