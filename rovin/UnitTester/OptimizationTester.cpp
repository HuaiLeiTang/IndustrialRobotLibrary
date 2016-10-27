#include <rovin\Optimizer\NonlinearOptimization.h>
#include <rovin\Math\Function.h>
#include <rovin\EnergyOptimization\PTPOptimization.h>
#include <rovin\Renderer\OSG_simpleRender.h>

#include "efortRobot.h"

#include <iostream>
#include <fstream>
#include <conio.h>
#include <time.h>
#include <string.h>

using namespace rovin;
using namespace std;

class TestObjFunction;
class TestEqFunction;
class TestIneqFunction;

void calculatePathfromeBSpline(MatrixX& result, const VectorX& t, BSpline<-1, -1, -1>& spline);
void calculatePathfromeBSpline(MatrixX& result, const unsigned int numOfdata, const Real ti, const Real tf, BSpline<-1, -1, -1>& spline);
void loadData_ys(MatrixX& data, std::string filename);
void saveMatrixX2txt(MatrixX in, std::string filename);
void saveVectorX2txt(VectorX in, std::string filename);
void calculateTorqueTrajectory(const MatrixX& q, const MatrixX& qdot, const MatrixX& qddot, MatrixX& torque);
Real makeRandLU(Real lower, Real upper)
{
	return lower + (upper - lower)*((double)rand() / 32767.0);
}

//std::string file = "D:/jkkim/Documents/matlabTest/opt/";
std::string file = "C:/Users/crazy/Desktop/Time optimization/nloptMMA test/";
std::string filePath = "C:/Users/crazy/Desktop/Time optimization/pathOptTest/";

SerialOpenChainPtr robot(new efortRobot());

int main()
{
	cout << "Test" << endl;
	cout << "Test Hong" << endl;

	unsigned int dof = robot->getNumOfJoint();

	StatePtr initState, finalState;
	initState = robot->makeState();
	finalState = robot->makeState();

	VectorX init_q(dof), init_qdot(dof), init_qddot(dof);
	VectorX final_q(dof), final_qdot(dof), final_qddot(dof);

	init_q << 1.0854, 1.02654, 0.798359, 2.97849, 1.50724, 1.45496;
	init_qdot.setZero(); init_qddot.setZero();
	final_q << -1.31465, 0.128787, -0.546992, 2.83671, -1.8583, 2.95029;
	final_qdot.setZero(); final_qddot.setZero();

	initState->setJointStatePos(init_q);initState->setJointStateVel(init_qdot);initState->setJointStateAcc(init_qddot);
	finalState->setJointStatePos(final_q);finalState->setJointStateVel(final_qdot);finalState->setJointStateAcc(final_qddot);

	///////////////////////////////// OPTIMIZATION /////////////////////////////////
	vector<bool> optJoint(robot->getNumOfJoint());
	optJoint[0] = optJoint[1] = optJoint[2] = true;
	Real tf = 5.0;
	int numOfOptCP = 4;
	int orderOfBSpline = 4;

	//MatrixX initM;
	//loadData_ys(initM, "C:/Users/crazy/Desktop/Time optimization/pathOptTest/dataOptResult2.txt");

	std::cout << "------- [NLOPT RESULT] -------" << endl;
	PTPOptimization PTPManagerNlopt(robot, optJoint, orderOfBSpline, numOfOptCP, 20, tf, initState, finalState, OptimizationType::nlopt , ObjectiveFunctionType::energyloss);
	//PTPManagerNlopt.initM = initM;
	PTPManagerNlopt.generateTrajectory();

	//BSpline<-1, -1, -1> qresult = PTPManagerNlopt._shared->_qSpline;
	//cout << qresult.getControlPoints() << endl;
	//cout << qresult.getKnots() << endl;
	//BSpline<-1, -1, -1> qdotresult = PTPManagerNlopt._shared->_qdotSpline;
	//BSpline<-1, -1, -1> qddotresult = PTPManagerNlopt._shared->_qddotSpline;
	//MatrixX qresultM, qdotresultM, qddotresultM;
	//calculatePathfromeBSpline(qresultM, 1000, 0, tf, qresult);
	//calculatePathfromeBSpline(qdotresultM, 1000, 0, tf, qdotresult);
	//calculatePathfromeBSpline(qddotresultM, 1000, 0, tf, qddotresult);
	//saveMatrixX2txt(qresultM, filePath + "datafromRovinq.txt");
	//saveMatrixX2txt(qdotresultM, filePath + "datafromRovinqdot.txt");
	//saveMatrixX2txt(qddotresultM, filePath + "datafromRovinqddot.txt");

	//std::cout << "------- [GCMMA RESULT] -------" << endl;
	//PTPOptimization PTPManagerManualOptTR(robot, optJoint, orderOfBSpline, numOfOptCP, 20, tf, initState, finalState, OptimizationType::GCMMA, ObjectiveFunctionType::effort);
	//PTPManagerManualOptTR.generateTrajectory();

	//std::cout << "------- [GCMMA TRUST_REGION RESULT] -------" << endl;
	//PTPOptimization PTPManagerManualOptTR(robot, optJoint, orderOfBSpline, numOfOptCP, 20, tf, initState, finalState, OptimizationType::GCMMA_TR, ObjectiveFunctionType::effort);
	//PTPManagerManualOptTR.generateTrajectory();

	//std::cout << "------- [GCMMA GD RESULT] -------" << endl;
	//PTPOptimization PTPManagerManualOptGD(robot, optJoint, orderOfBSpline, numOfOptCP, 20, tf, initState, finalState, OptimizationType::GCMMA_GD, ObjectiveFunctionType::effort);
	//PTPManagerManualOptGD.generateTrajectory();

	//MatrixX cpp(dof, 6 + numOfOptCP);
	//for (int i = 0; i < 3; i++)
	//{
	//	cpp.col(i) = PTPManagerManualOptTR._initialCP[i];
	//	cpp.col(6 + numOfOptCP - i) = PTPManagerManualOptTR._finalCP[i];
	//}
	//cpp.block(3, 3, 3, numOfOptCP) = PTPManagerManualOptTR._noptJointCP;
	//cout << cpp << endl;
	//for (int i = 0; i < 3; i++)
	//{
	//	for (int j = 0; j < numOfOptCP; j++)
	//	{
	//		cpp(i, 3 + j) = PTPManagerManualOptTR.initX(numOfOptCP * i + j);
	//	}
	//}
	//cout << cpp << endl;
	//cout << PTPManagerManualOptTR._knot << endl;

	///////////////////////////////// SAVE & RENDERING /////////////////////////////////
	//bool renderingswi = true;

	//PTPOptimization* PTPManager1 = &PTPManagerNlopt;
	//PTPOptimization* PTPManager2 = &PTPManagerManualOptTR;
	//int datanum = 2000;
	//Real stepsize = (tf - 0.0) / 2000;
	//Real t = 0.0;
	//BSpline<-1, -1, -1> nloptSpline = PTPManager1->_shared->_qSpline;
	//BSpline<-1, -1, -1> MMASpline = PTPManager2->_shared->_qSpline;
	//BSpline<-1, -1, -1> nloptqdotSpline = PTPManager1->_shared->_qdotSpline;
	//BSpline<-1, -1, -1> MMAqdotSpline = PTPManager2->_shared->_qdotSpline;
	//BSpline<-1, -1, -1> nloptqddotSpline = PTPManager1->_shared->_qddotSpline;
	//BSpline<-1, -1, -1> MMAqddotSpline = PTPManager2->_shared->_qddotSpline;
	//MatrixX nloptTraj(dof, datanum), MMATraj(dof, datanum);
	//MatrixX nloptqdotTraj(dof, datanum), MMAqdotTraj(dof, datanum);
	//MatrixX nloptqddotTraj(dof, datanum), MMAqddotTraj(dof, datanum);

	//for (int i = 0; i < datanum; i++)
	//{
	//	nloptTraj.col(i) = nloptSpline(t);
	//	MMATraj.col(i) = MMASpline(t);
	//	t += stepsize;
	//}
	//t = 0.0;
	//for (int i = 0; i < datanum; i++)
	//{
	//	nloptqdotTraj.col(i) = nloptqdotSpline(t);
	//	MMAqdotTraj.col(i) = MMAqdotSpline(t);
	//	t += stepsize;
	//}
	//t = 0.0;
	//for (int i = 0; i < datanum; i++)
	//{
	//	nloptqddotTraj.col(i) = nloptqddotSpline(t);
	//	MMAqddotTraj.col(i) = MMAqddotSpline(t);
	//	t += stepsize;
	//}

	//MatrixX nloptTorqueTraj(dof, datanum), MMATorqueTraj(dof, datanum);
	//calculateTorqueTrajectory(nloptTraj, nloptqdotTraj, nloptqddotTraj, nloptTorqueTraj);
	//calculateTorqueTrajectory(MMATraj, MMAqdotTraj, MMAqddotTraj, MMATorqueTraj);

	//saveMatrixX2txt(nloptTraj, file + "nlopt/q.txt");
	//saveMatrixX2txt(nloptqdotTraj, file + "nlopt/qdot.txt");
	//saveMatrixX2txt(nloptqddotTraj, file + "nlopt/qddot.txt");
	//saveMatrixX2txt(nloptTorqueTraj, file + "nlopt/torque.txt");

	//saveMatrixX2txt(MMATraj, file + "MMA/q.txt");
	//saveMatrixX2txt(MMAqdotTraj, file + "MMA/qdot.txt");
	//saveMatrixX2txt(MMAqddotTraj, file + "MMA/qddot.txt");
	//saveMatrixX2txt(MMATorqueTraj, file + "MMA/torque.txt");

	//// constraint save
	//VectorX qct(dof * 2), qdotct(dof * 2), qddotct(dof * 2), tct(dof * 2);
	//for (unsigned int i = 0; i < dof; i++)
	//{
	//	qct(i) = robot->getMotorJointPtr(i)->getLimitPosLower();
	//	qct(i + dof) = robot->getMotorJointPtr(i)->getLimitPosUpper();
	//	qdotct(i) = robot->getMotorJointPtr(i)->getLimitVelLower();
	//	qdotct(i + dof) = robot->getMotorJointPtr(i)->getLimitVelUpper();
	//	qddotct(i) = robot->getMotorJointPtr(i)->getLimitAccLower();
	//	qddotct(i + dof) = robot->getMotorJointPtr(i)->getLimitAccUpper();
	//	tct(i) = robot->getMotorJointPtr(i)->getLimitTorqueLower();
	//	tct(i + dof) = robot->getMotorJointPtr(i)->getLimitTorqueUpper();
	//}

	//saveVectorX2txt(qct, file + "qconstraint.txt");
	//saveVectorX2txt(qdotct, file + "qdotconstraint.txt");
	//saveVectorX2txt(qddotct, file + "qddotconstraint.txt");
	//saveVectorX2txt(tct, file + "torqueconstraint.txt");

	//if (renderingswi)
	//{
	//	OSG_simpleRender renderer(*robot, *initState, 800, 800);
	//	shared_ptr<Points> sp, ep;
	//	shared_ptr<Line> l_nlopt, l_MMA;
	//	Real psize = 20.0;

	//	SE3 T(Vector3(0.0, 0.0, 0.12));
	//	initState->setJointStatePos(init_q); robot->solveForwardKinematics(*initState);
	//	Vector3 spos = (initState->getLinkStateSE3(dof)*T).getPosition();
	//	finalState->setJointStatePos(final_q); robot->solveForwardKinematics(*finalState);
	//	Vector3 epos = (finalState->getLinkStateSE3(dof)*T).getPosition();
	//	sp = shared_ptr<Points>(new Points);
	//	ep = shared_ptr<Points>(new Points);
	//	sp->push_back(osg::Vec3(spos(0), spos(1), spos(2)));
	//	ep->push_back(osg::Vec3(epos(0), epos(1), epos(2)));
	//	sp->setSize(psize); ep->setSize(psize);
	//	sp->setColor(0.0f, 0.0f, 1.0f); ep->setColor(1.0f, 0.0f, 0.0f);
	//	renderer.addGeometry(*sp);
	//	renderer.addGeometry(*ep);

	//	Vector3 nloptPos, MMAPos;
	//	l_nlopt = shared_ptr<Line>(new Line);
	//	l_MMA = shared_ptr<Line>(new Line);
	//	for (int i = 0; i < datanum; i++)
	//	{
	//		initState->setJointStatePos(nloptTraj.col(i));
	//		robot->solveForwardKinematics(*initState);
	//		nloptPos = (initState->getLinkStateSE3(dof)*T).getPosition();
	//		l_nlopt->push_back(osg::Vec3(nloptPos(0), nloptPos(1), nloptPos(2)));
	//	}
	//	for (int i = 0; i < datanum; i++)
	//	{
	//		initState->setJointStatePos(MMATraj.col(i));
	//		robot->solveForwardKinematics(*initState);
	//		MMAPos = (initState->getLinkStateSE3(dof)*T).getPosition();
	//		l_MMA->push_back(osg::Vec3(MMAPos(0), MMAPos(1), MMAPos(2)));
	//	}
	//	l_nlopt->setColor(0.0f, 0.0f, 1.0f);
	//	l_MMA->setColor(1.0f, 0.0f, 0.0f);
	//	renderer.addGeometry(*l_nlopt);
	//	renderer.addGeometry(*l_MMA);

	//	Real startT = 0.0, endT = tf;
	//	Real frame = 500;
	//	Real dT = 1.0 / frame;

	//	renderer.getViewer().realize();
	//	int count = 0;
	//	Real pastT = clock();

	//	while (1)
	//	{
	//		if ((clock() - pastT) / 1e+3 >= dT)
	//		{
	//			pastT = clock();

	//			if (count >= nloptTraj.cols())
	//			{
	//				count = 0;
	//			}
	//			initState->setJointStatePos(nloptTraj.col(count));
	//			robot->solveForwardKinematics(*initState);
	//			count += 5;
	//		}
	//		renderer.updateFrame();
	//	}
	//}

	std::cout << "====== Program Complete ======" << endl;
	_getch();

	return 0;
}

void calculatePathfromeBSpline(MatrixX& result, const VectorX& t, BSpline<-1, -1, -1>& spline)
{
	int row = spline(0.1).rows(), col = t.size();
	result.resize(row, col);
	for (int i = 0; i < col; i++)
		result.col(i) = spline(t(i));
}
void calculatePathfromeBSpline(MatrixX& result, const unsigned int numOfdata, const Real ti, const Real tf, BSpline<-1, -1, -1>& spline)
{
	int row = spline(0.1).rows();
	result.resize(row, numOfdata);
	Real dt = (tf - ti) / (Real)(numOfdata - 1);
	Real t = ti;

	for (unsigned int i = 0; i < numOfdata; i++)
	{
		result.col(i) = spline(t);
		t += dt;
	}
}

void loadData_ys(MatrixX& data, std::string filename)
{
	ifstream input(filename);
	if (input.fail()) cout << "파일 열기 실패" << endl;
	else cout << "파일 열기 성공" << endl;

	Real tmp;
	unsigned int cnt = 0;
	while (!input.eof()) {
		input >> tmp;
		cnt++;
	}
	input.close();

	ifstream trajectory(filename);
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

void saveMatrixX2txt(MatrixX in, std::string filename)
{
	std::ofstream fout;
	fout.open(filename);

	for (int i = 0; i < in.cols(); i++)
	{
		for (int j = 0; j < in.rows(); j++)
			fout << in(j, i) << '\t';
		fout << std::endl;
	}

	fout.close();
}

void saveVectorX2txt(VectorX in, std::string filename)
{
	std::ofstream fout;
	fout.open(filename);

	for (int i = 0; i < in.size(); i++)
	{
		fout << in(i) << endl;
	}
	fout.close();
}

void calculateTorqueTrajectory(const MatrixX& q, const MatrixX& qdot, const MatrixX& qddot, MatrixX& torque)
{
	int dof = q.rows();
	int datanum = q.cols();
	StatePtr state = robot->makeState();

	for (int i = 0; i < datanum; i++)
	{
		state->setJointStatePos(q.col(i));
		state->setJointStateVel(qdot.col(i));
		state->setJointStateAcc(qddot.col(i));

		robot->solveInverseDynamics(*state);
		torque.col(i) = state->getJointStateTorque();
	}
}

class TestObjFunction : public Function
{
public:
	TestObjFunction() {}

	VectorX func(const VectorX& x) const
	{
		VectorX val(1);
		val(0) = x(0)*x(0) + x(1)*x(1);
		return val;
	}
};

class TestEqFunction : public Function
{
public:
	TestEqFunction() {}

	VectorX func(const VectorX& x) const
	{
		VectorX val(1);
		val(0) = x(0) + x(1) - 5;
		return val;
	}
};

class TestIneqFunction : public Function
{
public:
	TestIneqFunction() {}

	VectorX func(const VectorX& x) const
	{
		VectorX val(1);
		val(0) = -x(1) + 3;
		return val;
	}
};