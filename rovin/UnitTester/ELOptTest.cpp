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

Real makeRandLU(Real lower, Real upper)
{
	return lower + (upper - lower)*((double)rand() / 32767.0);
}

std::string filePath = "C:/Users/crazy/Desktop/Time optimization/ELOptTest/";
SerialOpenChainPtr robot(new efortRobot());
unsigned int dof = robot->getNumOfJoint();

int main()
{
	Real mass = 0;
	Inertia m(mass);
	SE3 TOOLTIP(Vector3(0, 0, 0.12));
	SE3 T = robot->_socLink[dof]._M;
	m.changeFrame(T* TOOLTIP); m += robot->_socLink[dof]._G; robot->_socLink[dof]._G = m;

	VectorX init_q(dof), init_qdot(dof), init_qddot(dof);
	VectorX final_q(dof), final_qdot(dof), final_qddot(dof);
	init_qdot.setZero(); init_qddot.setZero();
	final_qdot.setZero(); final_qddot.setZero();

	init_q << -0.995613, 0.21022, -0.287521, -2.25377, -2.17067, -4.97274;
	final_q << 1.80321, 1.31698, -0.4104, -2.31149, -2.06037, -4.99;

	StatePtr initState, finalState;
	initState = robot->makeState();
	finalState = robot->makeState();

	initState->setJointStatePos(init_q); initState->setJointStateVel(init_qdot); initState->setJointStateAcc(init_qddot);
	finalState->setJointStatePos(final_q); finalState->setJointStateVel(final_qdot); finalState->setJointStateAcc(final_qddot);

	vector<bool> optJoint(robot->getNumOfJoint());
	optJoint[0] = optJoint[1] = optJoint[2] = true;
	Real tf = 5.0;
	int numOfOptCP = 5;
	int orderOfBSpline = 4;
	int numOfSampling = 30;

	//MatrixX initcp(6, 6 + numOfOptCP);
	//VectorX knots(6 + numOfOptCP + orderOfBSpline);

	//initcp << -0.995613, -0.996761, -0.947534, -0.879466, -0.804541, -0.733156, -0.657312, -0.581468, -0.504065, -0.426505, -0.348531, -0.269955, -0.191352, -0.112535, -0.0332274,
	//	0.0463051, 0.128028, 0.21038, 0.293145, 0.37677, 0.461042, 0.544872, 0.629287, 0.714174, 0.798631, 0.883231, 0.966487, 1.04874, 1.12998, 1.20958,
	//	1.2881, 1.36544, 1.44156, 1.51566, 1.59064, 1.65962, 1.73524, 1.78671,
	//	1.80329, 1.80321,
	//	0.21022, 0.209839, 0.224902, 0.246285, 0.270273, 0.293722, 0.319363, 0.346034, 0.374508, 0.404697, 0.437138, 0.472498, 0.510732, 0.550395, 0.589491,
	//	0.626734, 0.662796, 0.697287, 0.730624, 0.763211, 0.79505, 0.82572, 0.855515, 0.88417, 0.911263, 0.937133, 0.961711, 0.985703, 1.00985, 1.03486,
	//	1.06181, 1.09108, 1.12299, 1.15695, 1.19457, 1.2317, 1.2752, 1.30648, 1.31701, 1.31698,
	//	-0.287521, -0.287461, -0.290532, -0.294564, -0.298829, -0.302664, Q-0.306458, -0.309855, -0.312841, -0.315199, -0.316772, -0.317345, -0.316852, -0.315899, -0.315233,
	//	-0.315162, -0.315614, -0.316286, -0.316821, -0.317058, -0.316932, -0.316494, -0.315858, -0.315227, -0.314854, -0.314981, -0.315846, -0.317661, -0.320651, -0.324985,
	//	-0.330868, -0.338297, -0.347274, -0.357478, -0.369373, -0.381488, -0.396063, -0.406756, -0.410408, -0.4104,
	//	-2.25377, -2.25375, -2.25451, -2.25558, -2.25679, -2.25798, -2.25929, -2.26065, -2.26211, -2.26368, -2.26537, -2.26723, -2.26925, -2.27135, -2.27341,
	//	-2.27536, -2.27723, -2.27902, -2.28073, -2.28241, -2.28407, -2.28569, -2.28729, -2.28889, -2.29047, -2.29204, -2.29359, -2.29513, -2.29668, -2.29823,
	//	-2.2998, -2.3014, -2.30303, -2.30467, -2.30638, -2.308, -2.30981, -2.31108,
	//	-2.31149, -2.31149,
	//	-2.17067, -2.17071, -2.16926, -2.16721, -2.16489, -2.16262, -2.16013, -2.15752, -2.15472, -2.15174, -2.1485, -2.14495, -2.14109, -2.13707, -2.13314,
	//	-2.12941, -2.12583, -2.12243, -2.11915, -2.11593, -2.11277, -2.10968, -2.10661, -2.10356, -2.10055, -2.09754, -2.09458, -2.09163, -2.08867, -2.08571,
	//	-2.0827, -2.07964, -2.07653, -2.0734, -2.07013, -2.06704, -2.06357, -2.06116, -2.06037, -2.06037,
	//	-4.97274, -4.97273, -4.97296, -4.97328, -4.97364, -4.974, -4.97439, -4.9748,
	//	-4.97524, -4.9757, -4.97621, -4.97676, -4.97737, -4.978, -4.97861, -4.9792,
	//	-4.97976, -4.98029, -4.9808, -4.98131, -4.9818, -4.98228, -4.98276, -4.98324, -4.98371, -4.98418, -4.98465, -4.98511, -4.98557, -4.98603, -4.98651,
	//	-4.98698, -4.98747, -4.98796, -4.98847, -4.98896, -4.9895, -4.98988, -4.99,
	//	-4.99;
	//knots << 0, 0, 0, 0, 0.135135, 0.27027, 0.405405, 0.540541, 0.675676, 0.810811, 0.945946, 1.08108, 1.21622, 1.35135, 1.48649, 1.62162, 1.75676,
	//	1.89189, 2.02703, 2.16216, 2.2973, 2.43243, 2.56757, 2.7027, 2.83784, 2.97297, 3.10811, 3.24324, 3.37838, 3.51351, 3.64865, 3.78378, 3.91892,
	//	4.05405, 4.18919, 4.32432, 4.45946, 4.59459, 4.72973, 4.86486, 5,
	//	5, 5, 5;
	//BSpline<-1, -1, -1> initBSpline(knots, initcp);

	PTPOptimization PTPManagerNlopt(robot, optJoint, orderOfBSpline, numOfOptCP, numOfSampling, tf, initState, finalState, OptimizationType::nlopt, ObjectiveFunctionType::energyloss);
	//PTPOptimization PTPManagerNlopt(robot, optJoint, initBSpline, numOfSampling, initState, finalState, OptimizationType::nlopt, ObjectiveFunctionType::energyloss);
	PTPManagerNlopt.generateTrajectory();


	cout << "Program Complete" << endl;
	_getch();
	return 0;
}