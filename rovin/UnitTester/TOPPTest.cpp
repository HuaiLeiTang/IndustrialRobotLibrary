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
	TOPP topp(q_data, robot, ds, vi, vf, si, sf);
	topp.generateTrajectory();
	//topp.saveMVCandSP2txt();
	//backwardIntegrationTest(topp);

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	//std::vector<Real> min;
	//std::vector<Real> max;
	//VectorX qs(topp._dof);
	//VectorX vec(topp._dof * 2);
	//VectorX qconts = topp._velConstraint;
	//qconts *= 4;
	//VectorX left_vec(topp._dof * 2);

	//Vector2 result;

	//for (int i = 0; i < topp._dof; i++)
	//{
	//	cout << robot->getMotorJointPtr(i)->getLimitVelUpper() << endl;
	//	//cout << "[ joint " << "i" << " qdot max min ]" << endl;
	//	//cout << "max : " << robot->getMotorJointPtr(i)->getLimitVelUpper() << endl;
	//	//cout << "min : " << robot->getMotorJointPtr(i)->getLimitVelLower() << endl;
	//}

	//cout << topp._dqds(0.1) << endl;

	//Real s_cur = si;
	//while (s_cur < sf)
	//{
	//	qs = topp._dqds(s_cur);
	//	
	//	for (int i = 0; i < topp._dof; i++)
	//	{
	//		vec(i) = qs(i);
	//		vec(i + topp._dof) = -qs(i);
	//		left_vec(i) = qconts(i);
	//		left_vec(i + topp._dof) = -qconts(i + topp._dof);
	//	}

	//	result[0] = -std::numeric_limits<Real>::max();
	//	result[1] = std::numeric_limits<Real>::max();
	//	for (int i = 0; i < topp._dof * 2; i++)
	//	{
	//		Real tmp = left_vec(i) / vec(i);
	//		if (vec[i] > RealEps) // upper bound beta
	//		{
	//			if (tmp < result[1])
	//				result[1] = tmp;
	//		}
	//		else if (vec[i] < -RealEps)// lower bound alpha
	//		{
	//			if (tmp > result[0])
	//				result[0] = tmp;
	//		}
	//	}
	//	result[0] = std::max(result[0], 0.0);

	//	min.push_back(result[0]);
	//	max.push_back(result[1]);

	//	s_cur += topp._ds;
	//}
	//topp.saveRealVector2txt(min, "C:/Users/crazy/Desktop/Time optimization/qdot minmax result/min.txt");
	//topp.saveRealVector2txt(max, "C:/Users/crazy/Desktop/Time optimization/qdot minmax result/max.txt");
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////


	cout << "Program complete" << endl;
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

void backwardIntegrationTest(TOPP& topp)
{
	topp.calcMVC();
	int cnt = 0;
	for (int i = 0; i < topp.s_MVC_jk.size(); i = i + 5)
	{
		if (topp.s_MVC_jk[i] > 0.400)
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