/*!
 *	\file	PTPOptimizer.h
 *	\date	2016.02.18
 *	\author	Youngsuk (crazyhys@gmail.com)
 *	\brief	PTPOptimization class
*/

#pragma once

#include <rovin\Optimizer\NonlinearOptimization.h>
#include <rovin\Dynamics\SerialOpenChain.h>
#include <rovin\Math\GaussianQuadrature.h>
#include <rovin\Math\Interpolation.h>

// #ifndef
// #define EPS 100
// #endif
// 이런식으로 쓰면 될까???

namespace rovin {

	class PTPOptimizer;
	class sharedResource;

	class effortFunction;
	class NonlinearInequalityConstraint;

	typedef std::shared_ptr<sharedResource> sharedResourcePtr;



	class PTPoptimizer
	{
		friend class sharedResource;
		friend class effortFunction;
		friend class NonlinearInequalityConstraint;
	public:
		enum ObjectiveFunctionType {effort, energyloss};
	
	private:
		NonlinearOptimization _optimizer;
		std::vector<sharedResourcePtr> _shared;

		// functions
		FunctionPtr _objectFunc;
		FunctionPtr _IneqFunc;
		FunctionPtr _linearIneqFunc;
		FunctionPtr _nonlinearIneqFunc;

		SerialOpenChainPtr _soc;

		std::vector<bool> _optJoint;
		std::vector<unsigned int> _optJointIdx;
		std::vector<unsigned int> _noptJointIdx;
		unsigned int _numOfOptJoint;

		unsigned int _orderOfBSpline;
		unsigned int _numOfOptCP;
		unsigned int _numOfGQSample;

		//unsigned int _numOfWayPoint;
		unsigned int _numOfShared;

		std::vector<Real> _tf;

		std::vector<StatePtr> _constraintState;

		// Gaussian quadrature
		GaussianQuadrature GQ;

	public:
		PTPoptimizer(const SerialOpenChainPtr& soc, const std::vector<bool>& optJoint, 
			const unsigned int orderOfBSpline, const unsigned int numOfOptCP, const unsigned int numOfGQSample, 
			const std::vector<Real>& tf, const std::vector<StatePtr>& constraintState);

		// Optimization function
		void makeObjectiveFunction();
		void makeIneqConstraintFunction();

		// calculate energy efficiency trajectory
		void generateTrajectory();
	};

	class sharedResource
	{
	private:

		PTPOptimizer* _PTPOptimizer;

		VectorX _params;

		std::vector<VectorX> _tau;
		std::vector<MatrixX> _dtaudp;
		std::vector<StatePtr> _state;

		std::vector<MatrixX> _dqdp;
		std::vector<MatrixX> _dqdotdp;
		std::vector<MatrixX> _dqddotdp;

	public:
		BSpline<-1, -1, -1> _qSpline;
		BSpline<-1, -1, -1> _qdotSpline;
		BSpline<-1, -1, -1> _qddotSpline;

		MatrixX _dPdP;
		MatrixX _dQdP;
		MatrixX _dRdP;

		std::vector<VectorX> _P;
		std::vector<VectorX> _Q;
		std::vector<VectorX> _R;

		// B spline variables
		VectorX _knot;
		std::vector<VectorX> _initialCP;
		std::vector<VectorX> _finalCP;
		MatrixX _noptJointCP;

	public:
		sharedResource(PTPOptimizer* PTPOptimizer);

		void makeBSpline(const VectorX& params);


	};

}



