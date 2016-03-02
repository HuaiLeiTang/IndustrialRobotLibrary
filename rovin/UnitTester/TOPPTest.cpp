#include <iostream>
#include <fstream>
#include <conio.h>
#include <rovin\TimeOptimization\TOPP.h>

#include "efortRobot.h"
#include <string>

using namespace std;
using namespace rovin;

void loadData(MatrixX& data);
void backwardIntegrationTest(TOPP& topp);
void savessdotResult(TOPP& topp);
void saveTorqueResult(TOPP& topp);
void findswitchingPoint(TOPP& topp);
void savevelminmax(TOPP& topp);

SerialOpenChainPtr robot(new efortRobot());
StatePtr state;
unsigned int dof;

int main()
{
	state = robot->makeState();
	dof = robot->getNumOfJoint();

	MatrixX q_data;
	loadData(q_data);

	// Time optimization
	Real ds = 1e-3, vi = 0, vf = 0, si = 0, sf = 1;
	TOPP topp(q_data, robot, ds, vi, vf, si, sf, CONSTRAINT_TYPE::TORQUE_ACC);
	//topp.saveTrap2txt();
	//topp.generateTrajectory();
	//savessdotResult(topp);
	//topp.saveMVCandSP2txt();
	//saveTorqueResult(topp);
	//cout << "a" << endl;
	//findswitchingPoint(topp);

	//topp.generateTrajectory();
	//savessdotResult(topp);
	//topp.saveMVCandSP2txt();
	//saveTorqueResult(topp);


	cout << "Final time : " << topp.getFinalTime() << endl;
	cout << "Program complete" << endl;
	_getch();
	return 0;
}

void savevelminmax(TOPP& topp) //--> minmax 함수 바꾸기
{
	std::vector<Real> min;

	Real s_cur = 0, sdot_cur;
	VectorX qs(topp._dof);
	VectorX qmax(topp._dof);
	for (int i = 0; i < topp._dof; i++)
	{
		qmax(i) = topp._velConstraint(i);
		qmax(i) *= qmax(i);
	}
	while (s_cur < topp._sf)
	{
		qs = topp._dqds(s_cur);
		for (int i = 0; i < topp._dof; i++)
			qs(i) *= qs(i);

		Real tmp_max = std::numeric_limits<Real>::max();
		for (int i = 0; i < topp._dof; i++)
		{
			if (qmax(i) / qs(i) < tmp_max)
				tmp_max = qmax(i) / qs(i);
		}
		sdot_cur = sqrt(tmp_max);
		min.push_back(sdot_cur);
		s_cur += topp._ds;
	}
	topp.saveRealVector2txt(min, "C:/Users/crazy/Desktop/Time optimization/min.txt");
}

void loadData(MatrixX& data)
{

	//ifstream input("trajectory.txt");
	//ifstream input("C:/Users/crazy/Desktop/Time optimization/trajectory text/trajectory_wy.txt");
	//ifstream input("C:/Users/ksh/Documents/trajectory_wy.txt");
	//ifstream input("D:/jkkim/Documents/trajectory_wy.txt");
	ifstream input("C:/Users/crazy/Desktop/Time optimization/trajectory text/trajectory_wy.txt");
	//ifstream input("C:/Users/ksh/Documents/trajectory_wy.txt");

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
	//ifstream trajectory("C:/Users/ksh/Documents/trajectory_wy.txt");
	//ifstream trajectory("D:/jkkim/Documents/trajectory_wy.txt");
	ifstream trajectory("C:/Users/crazy/Desktop/Time optimization/trajectory text/trajectory_wy.txt");
	//ifstream trajectory("C:/Users/ksh/Documents/trajectory_wy.txt");

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

void backwardIntegrationTest(TOPP& topp)
{
	topp.calcMVC();
	int cnt = 0;
	for (int i = 0; i < topp.s_MVC_jk.size(); i = i + 5)
	{
		//if (topp.s_MVC_jk[i] > 0.400)
		if (topp.s_MVC_jk[i] > 0.370)
		{
			string s_st = "C:/Users/crazy/Desktop/Time optimization/tmp/s";
			string sd_st = "C:/Users/crazy/Desktop/Time optimization/tmp/sd";
			Real s_cur = topp.s_MVC_jk[i];
			Real sdot_cur = topp.sd_MVC_jk[i];
			Vector2 alphabeta;
			std::vector<Real> s_result;
			std::vector<Real> sd_result;
			for (int j = 0; j < 400; j++)
			{
				s_result.push_back(s_cur);
				sd_result.push_back(sdot_cur);
				alphabeta = topp.determineAlphaBeta(s_cur, sdot_cur);
				topp.backwardIntegrate(s_cur, sdot_cur, alphabeta(0));
			}
			//cout << "cnt : " << cnt << endl;
			s_st += to_string(cnt);
			s_st += ".txt";
			sd_st += to_string(cnt);
			sd_st += ".txt";
			topp.saveRealVector2txt(s_result, s_st);
			topp.saveRealVector2txt(sd_result, sd_st);
			s_result.clear();
			sd_result.clear();
			cnt++;
		}
	}
}

void savessdotResult(TOPP& topp)
{
	std::vector<Real> s_result;
	std::vector<Real> sdot_result;
	std::list<Real>::const_iterator s_it = topp.gets().begin();
	std::list<Real>::const_iterator sdot_it = topp.getsdot().begin();
	for (s_it = topp.gets().begin(); s_it != topp.gets().end(); ++s_it)
	{
		s_result.push_back(*(s_it));
		sdot_result.push_back(*(sdot_it));
		sdot_it++;
	}
	topp.saveRealVector2txt(s_result, "C:/Users/crazy/Desktop/Time optimization/s_result.txt");
	topp.saveRealVector2txt(sdot_result, "C:/Users/crazy/Desktop/Time optimization/sdot_result.txt");
	//saveRealVector2txt(s_result, "D:/jkkim/Documents/matlabTest/s_result.txt");
	//saveRealVector2txt(sdot_result, "D:/jkkim/Documents/matlabTest/sdot_result.txt");
	//topp.saveRealVector2txt(s_result, "C:/Users/ksh/Documents/MATLAB/s_result.txt");
	//topp.saveRealVector2txt(sdot_result, "C:/Users/ksh/Documents/MATLAB/sdot_result.txt");
}

void saveTorqueResult(TOPP& topp)
{
	MatrixX torque = topp.getTorqueTrajectory();
	std::vector<vector<Real>> torque_vec(topp.getdof());

	for (int i = 0; i < topp.getdof(); i++)
	{
		string torque_st = "C:/Users/crazy/Desktop/Time optimization/torque";
		for (int j = 0; j < torque.row(0).size(); j++)
			torque_vec[i].push_back(torque(i, j));

		torque_st += to_string(i+1);
		torque_st += ".txt";
		topp.saveRealVector2txt(torque_vec[i], torque_st);
	}
}

void findswitchingPoint(TOPP& topp)
{
	Real s_cur = 0;
	bool swi = false;
	while (s_cur < topp._sf)
	{
		swi = topp.findNearestSwitchPoint(s_cur);
		s_cur = topp._switchPoint[topp._switchPoint.size() - 1]._s;
		if (!swi)
			break;
	}

	for (int i = 0; i < topp._switchPoint.size(); i++)
	{
		cout << i << "-th switching point" << endl;
		cout << "s : " << topp._switchPoint[i]._s << endl;
		cout << "sdot : " << topp._switchPoint[i]._sdot << endl;
		cout << "id : " << topp._switchPoint[i]._id << endl;
 	}
}