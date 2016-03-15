#include <iostream>
#include <fstream>
#include <conio.h>
#include <rovin\TimeOptimization\TOPP.h>
#include <rovin\TimeOptimization\P2PTimeOptTemp.h>
//#include <rovin\Renderer\OSG_simpleRender.h>

#include "efortRobot.h"
#include <string>
#include <time.h>

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
	//MatrixX q_data;
	//loadData(q_data);

	//Real ds = 1e-3, vi = 0, vf = 0, si = 0, sf = 1;

	//TOPP topp(q_data, robot, vi, vf, ds, si, sf, CONSTRAINT_TYPE::TORQUE_VEL_ACC);
	//topp.generateTrajectory();
	//cout << "[ consider only torque constraint ]" << endl;
	//cout << "Final time : " << topp.getFinalTime() << endl << endl;

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
	//VectorX qvel;
	//bool forward;

	//Vertex* nVertex = new Vertex();
	//nVertex->setconfig(qnear);
	//nVertex->setconfigVel(qnearVel);

	//AVP_RRT avp_rrt(robot, CONSTRAINT_TYPE::TORQUE_VEL_ACC);
	//avp_rrt.interpolate(nVertex, qrand, dist, forward, Pnew, qnew, qvel);

	//cout << "[qnew]" << endl;
	//cout << qnew << endl << endl;

	//delete nVertex;

	////////////////////////////////////////////////////////////////////////
	// AVP test
	state = robot->makeState();
	dof = robot->getNumOfJoint();

	AVP_RRT avp_rrt(robot, CONSTRAINT_TYPE::TORQUE);

	VectorX qs(dof), qg(dof), qsdot(dof), qgdot(dof);
	qs << 0.115499, 0.220374, 0.151515, 0.15633, 0.0438677, 0.0182414;
	qg << 1.15739, 1.11361, 1.04554, 1.25135, 1.0009, 0.790497;
	qsdot.setZero(); qgdot.setZero();
	
	VectorX q1(dof), q1dot(dof);
	q1 << 0.638163, 0.604759, 0.800566, 0.846903, 0.488474, 0.422513;
	q1dot.setZero();

	// experiment joint angle
	VectorX qex1(dof), qex2(dof);
	qex1 << 0.359377, 0.314241, 0.139612, 0.258737, -0.0754488, -0.377798;
	qex2 << 0.425051, 0.250989, 0.244523, -0.188007, -0.0782789, 0.123605;

	std::vector<WayPoint> wayPt(2);
	wayPt[0].setJointq(qs); wayPt[0].setJointqdot(qsdot);
	wayPt[1].setJointq(q1); wayPt[1].setJointqdot(q1dot);

	avp_rrt.setWayPoints(wayPt);
	avp_rrt.treeInitialization(0);
	cout << "start tree vertex config : " << avp_rrt._startTree._nodes[0]->_config << endl;
	cout << "start tree vertex configvel : " << avp_rrt._startTree._nodes[0]->_configVel << endl;
	cout << "end tree vertex config : " << avp_rrt._goalTree._nodes[0]->_config << endl;
	cout << "end tree vertex configvel : " << avp_rrt._goalTree._nodes[0]->_configVel << endl;

	Vertex * nVertex = avp_rrt._startTree._nodes[0];
	//VectorX qrand = q1;
	VectorX qrand = qg;
	//VectorX qrand = qex2;
	Real dist = 1.0;;
	bool forward = true;
	std::list<VectorX> Pnew; 
	VectorX qnew;
	VectorX qvel;

	avp_rrt.interpolate_tmp(nVertex, qrand, dist, forward, Pnew, qnew, qvel);

	cout << "Pnew initial : " << Pnew.front() << endl;
	cout << "Pnew final : " << Pnew.back() << endl;

	MatrixX q_data(Pnew.front().size(), Pnew.size());
	unsigned int cnt = 0;
	for (std::list<VectorX>::iterator it = Pnew.begin(); it != Pnew.end(); it++)
		q_data.col(cnt++) = *it;
	int data_num = Pnew.size();

	//OSG_simpleRender renderer(*robot, *state, 600, 600);
	//Real spsize = 50.0;

	//std::shared_ptr<Points> startpoint, endpoint;
	//std::shared_ptr<Line> line;

	//state->setJointStatePos(qs);
	//robot->solveForwardKinematics(*state);
	//Vector3 startpos = robot->calculateEndeffectorFrame(state).getPosition();

	//state->setJointStatePos(qrand);
	//robot->solveForwardKinematics(*state);
	//Vector3 endpos = robot->calculateEndeffectorFrame(state).getPosition();
	//
	//startpoint = shared_ptr<Points>(new Points);
	//endpoint = shared_ptr<Points>(new Points);
	//startpoint->push_back(osg::Vec3(startpos(0), startpos(1), startpos(2)));
	//endpoint->push_back(osg::Vec3(endpos(0), endpos(1), endpos(2)));
	//startpoint->setSize(spsize);
	//endpoint->setSize(spsize);
	//startpoint->setColor(0.0f, 0.0f, 1.0f);
	//endpoint->setColor(1.0f, 0.0f, 0.0f);
	//renderer.addGeometry(*startpoint);
	//renderer.addGeometry(*endpoint);

	//Vector3 robotpos;
	//line = shared_ptr<Line>(new Line);
	//for (unsigned int i = 0; i < q_data.cols(); i++)
	//{
	//	state->setJointStatePos(q_data.col(i));
	//	robot->solveForwardKinematics(*state);
	//	robotpos = robot->calculateEndeffectorFrame(state).getPosition();
	//	line->push_back(osg::Vec3(robotpos(0), robotpos(1), robotpos(2)));
	//}
	//renderer.addGeometry(*line);

	//renderer.getViewer().realize();
	//double frameRate = 10;

	//cnt = 0;
	//double c = clock();
	//while (1)
	//{
	//   if (clock() - c >= 100 / frameRate)
	//   {
	//      if (cnt == data_num) cnt = 0;
	//	  state->setJointStatePos(q_data.col(cnt));
	//      robot->solveForwardKinematics(*state);
	//	  cnt++;
	//   }
	//   renderer.updateFrame();
	//}

	std::cout << "[ TOPP 수행 ]" << std::endl;
	Real ds = 0.5e-2, vi = 0, vf = 0, si = 0, sf = 1;
	TOPP topp(q_data, robot, vi, vf, ds, si, sf, CONSTRAINT_TYPE::TORQUE);

	std::vector<Vector2, Eigen::aligned_allocator<Vector2>> allMVCPoints;
	std::vector<SwitchPoint> allSwitchPoints;
	topp.calculateAllMVCPoint();
	topp.calculateAllSwitchPoint();
	allMVCPoints = topp.getAllMVCPoint();
	allSwitchPoints = topp.getAllSwitchPoint();

	cout << "number of switch points : " << allSwitchPoints.size() << endl;

	topp.saveMVC(allMVCPoints);
	topp.saveSwitchPoint(allSwitchPoints);

	topp.generateTrajectory();
	topp.saveRealList2txt(topp._s, "C:/Users/crazy/Desktop/Time optimization/avp test/s_result.txt");
	topp.saveRealList2txt(topp._sdot, "C:/Users/crazy/Desktop/Time optimization/avp test/sdot_result.txt");


	cout << "Program complete" << endl;
	_getch();
	return 0;
}

void loadData(MatrixX& data)
{
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