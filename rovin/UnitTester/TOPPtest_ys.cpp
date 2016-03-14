#include <iostream>
#include <fstream>
#include <conio.h>
#include <rovin\TimeOptimization\TOPP.h>
#include <rovin\TimeOptimization\P2PTimeOptTemp.h>
#include <rovin\Renderer\OSG_simpleRender.h>

#include "efortRobot.h"
#include <string>
#include <time.h>

using namespace std;
using namespace rovin;

void loadData(MatrixX& data);

std::shared_ptr<Points> startpoint, endpoint;

SerialOpenChainPtr robot(new efortRobot());
StatePtr state;
unsigned int dof;

int main()
{
	////////////////////////////////////////////////////////////////////////
	// TOPP test
	MatrixX q_data;
	loadData(q_data);

	Real ds = 1e-3, vi = 0, vf = 0, si = 0, sf = 1;

	TOPP topp(q_data, robot, vi, vf, ds, si, sf, CONSTRAINT_TYPE::TORQUE_VEL_ACC);
	topp.generateTrajectory();
	cout << "[ consider only torque constraint ]" << endl;
	cout << "Final time : " << topp.getFinalTime() << endl << endl;

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
	//MatrixX q_data;
	//loadData(q_data);
	//VectorX qnear = q_data.col(0);
	//VectorX qnearVel = VectorX::Ones(6) * 0.1;
	//VectorX qrand = q_data.col(30);
	//Real dist = 1.0;

	//cout << "[qnear]" << endl;
	//cout << qnear << endl << endl;
	//cout << "[qrand]" << endl;
	//cout << qrand << endl << endl;
	//cout << "[qnearVel]" << endl;
	//cout << qnearVel << endl << endl;

	//std::list<VectorX> Pnew;
	//VectorX qnew;

	//Vertex* nVertex = new Vertex();
	//nVertex->setconfig(qnear);
	//nVertex->setconfigVel(qnearVel);

	//AVP_RRT avp_rrt(robot, CONSTRAINT_TYPE::TORQUE_VEL_ACC);
	//avp_rrt.interpolate(nVertex, qrand, dist, Pnew, qnew);

	//cout << "[qnew]" << endl;
	//cout << qnew << endl << endl;

	//delete nVertex;

	////////////////////////////////////////////////////////////////////////
	// rendering
	//state = robot->makeState();
	//MatrixX q_data;
	//loadData(q_data);
	//int data_num = q_data.cols();

	//state->setJointStatePos(0, q_data(0, 0));state->setJointStatePos(1, q_data(1, 0));
	//state->setJointStatePos(2, q_data(2, 0));state->setJointStatePos(3, q_data(3, 0));
	//state->setJointStatePos(4, q_data(4, 0));state->setJointStatePos(5, q_data(5, 0));
	//robot->solveForwardKinematics(*state);

	//OSG_simpleRender renderer(*robot, *state, 600, 600);
	//renderer.getViewer().realize();

	//Real spsize = 100.0;
	//startpoint = shared_ptr<Points>(new Points);
	//Vector3 startpos = robot->calculateEndeffectorFrame(state).getPosition();
	////startpoint->push_back(osg::Vec3(startpos(0), startpos(1), startpos(2)));
	//startpoint->push_back(osg::Vec3(0.0, 0.0, 0.0));
	//startpoint->setSize(spsize);
	//startpoint->setColor(1.0f, 1.0f, 01.0f);
	//renderer.addGeometry(*startpoint);

	//double frameRate = 30;

	//int cnt = 0;
	//double c = clock();
	//while (1)
	//{
	//   if (clock() - c >= 1000 / frameRate)
	//   {
	//      if (cnt == data_num) cnt = 0;
	//      state->setJointStatePos(0, q_data(0, cnt));
	//      state->setJointStatePos(1, q_data(1, cnt));
	//      state->setJointStatePos(2, q_data(2, cnt));
	//      state->setJointStatePos(3, q_data(3, cnt));
	//      state->setJointStatePos(4, q_data(4, cnt));
	//      state->setJointStatePos(5, q_data(5, cnt));
	//      cnt++;
	//      robot->solveForwardKinematics(*state);
	//   }
	//   renderer.updateFrame();
	//}

	cout << "Program complete" << endl;
	_getch();
	return 0;
}

void loadData(MatrixX& data)
{
	//ifstream input("D:/jkkim/Documents/trajectory_wy.txt");
	ifstream input("C:/Users/crazy/Desktop/Time optimization/trajectory text/trajectory_wy.txt");
	if (input.fail()) cout << "���� ���� ����" << endl;
	else cout << "���� ���� ����" << endl;

	Real tmp;
	unsigned int cnt = 0;
	while (!input.eof()) {
		input >> tmp;
		cnt++;
	}
	input.close();

	//ifstream trajectory("D:/jkkim/Documents/trajectory_wy.txt");
	ifstream trajectory("C:/Users/crazy/Desktop/Time optimization/trajectory text/trajectory_wy.txt");
	if (trajectory.fail()) cout << "���� ���� ����" << endl;
	else cout << "���� ���� ����" << endl;

	unsigned int dof = 6;
	unsigned int data_num = cnt / dof;
	MatrixX q(dof, data_num);

	for (unsigned int i = 0; i < data_num; i++)
		for (unsigned int j = 0; j < dof; j++)
			trajectory >> q(j, i);

	trajectory.close();

	data = q;
}