/*!
*	\file	PTPWayPointOptimization.h
*	\date	2016.02.17
*	\author	Youngsuk (crazyhys@gmail.com)
*	\brief	PTPWayPointOptimization class
*          this class make a energy efficiency optimization path of PToP problem with way points
*/

#pragma once

#include <rovin/Optimizer/NonlinearOptimization.h>

#include <rovin/Dynamics/SerialOpenChain.h>
#include <rovin/Math/GaussianQuadrature.h>
#include <rovin/Math/Interpolation.h>

namespace rovin {

	class PTPWayPointOptimizer;

	class PTPWayPointOptimization;
	class sharedResourceWayPoint;


	class effortFunctionWayPoint;
	class NonlinearInequalityConstraintWayPoint;

	//class energyLossFunction;
	//class EqualityConstraint;
	//class InequalityConstraint;

	typedef std::shared_ptr<sharedResourceWayPoint> sharedResourceWayPointPtr;
	typedef std::shared_ptr<PTPWayPointOptimization> PTPWayPointOptimizationPtr;
	

	class PTPWayPointOptimizer
	{
	private:
		PTPWayPointOptimizationPtr _wayPointOptimizer;
		
		unsigned int _numOfBspline;

		// result variable : b-spline
		std::vector<BSpline<-1, -1, -1>> _bsplineResult;

	public:
		PTPWayPointOptimizer(const SerialOpenChainPtr& soc, const std::vector<bool>& optJoint, const unsigned int orderOfBSpline, const unsigned int numOfCP, const unsigned int numOfGQSample, const std::vector<Real>& tf, const std::vector<StatePtr>& constraintState);
		void generateTrajectory();


		
	};



	class PTPWayPointOptimization
	{
		friend class sharedResourceWayPoint;
		friend class effortFunctionWayPoint;
		friend class NonlinearInequalityConstraintWayPoint;

	public:
		enum ObjectiveFunctionType { effort, energyloss };
		enum ConstraintCondition
		{
			INITIAL_ALL_FINAL_ALL, ///< initial pos, vel, acc constraint, final pos, vel, acc constraint
			INITIAL_ALL_FINAL_POS ///< initial pos, vel, acc constraint, final pos acc constraint
		};

	public:
		NonlinearOptimization _optimizer;
		sharedResourceWayPointPtr _shared;

		ConstraintCondition _constraintCondition;

		// functions
		FunctionPtr _objectFunc;
		FunctionPtr _IneqFunc;
		FunctionPtr _linearIneqFunc;
		FunctionPtr _nonlinearIneqFunc;
		//equalityConstraint _equlfunc;

		SerialOpenChainPtr _soc; ///< Serial open chain robot

		std::vector<bool> _optJoint;
		unsigned int _numOfOptJoint;
		std::vector<unsigned int> _optJointIdx;
		std::vector<unsigned int> _noptJointIdx;

		unsigned int _numOfConstraintCP;

		unsigned int _orderOfBSpline; ///< order of B spline
		unsigned int _numOfOptCP; ///< number of B spline control points for optimization except contraint control points
		unsigned int _numOfGQSample; ///< number of Gaussian quadrature sampling

		Real _tf; ///< traveling time

		StatePtr _initialState; ///< initial joint position
		StatePtr _finalState; ///< final joint position

		// B spline variables
		VectorX _knot; ///< number of knots = number of control points + order + 1
		std::vector<VectorX> _initialCP;
		std::vector<VectorX> _finalCP;
		MatrixX _noptJointCP;

		// Gaussian quadrature
		GaussianQuadrature GQ;

	public:
		PTPWayPointOptimization(const SerialOpenChainPtr& soc, const std::vector<bool>& optJoint, const unsigned int orderOfBSpline,
			const unsigned int numOfCP, const unsigned int numOfGQSample, const Real tf, 
			const StatePtr& initialState, const StatePtr& finalState, ConstraintCondition constraintCondition);

		// set function
		void setInitialState(const StatePtr& initState);
		void setFinalState(const StatePtr& finalState);
		void setConstraintCondition(const ConstraintCondition constraintCondition);
		void setFinalTime(const Real tf);


		// get function
		// Bspline 어케 return? BSpline<-1, -1, -1>&?
		// Bspline 결과 augment 시킬 수 있나?
		const BSpline<-1, -1, -1>& getResultBSpline() const;

		//
		void makeNonOptJointCP();

		// B spline function
		void makeBSplineKnot();
		void makeBoundaryCondition();

		// Optimization function
		void makeObjectiveFunction();
		void makeIneqConstraintFunction();

		// calculate energy efficiency trajectory
		void generateTrajectory();
	};

	// calculate inverse 
	class sharedResourceWayPoint
	{

	private:

		PTPWayPointOptimization* _PTPOptimizer;

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

	public:
		sharedResourceWayPoint(PTPWayPointOptimization* PTPOptimizer);

		void makeBSpline(const VectorX& params);

		void update(const VectorX& params);

		const std::vector<VectorX>& gettau(const VectorX& params);
		const std::vector<MatrixX>& getdtaudp(const VectorX& params);

	};

	class effortFunctionWayPoint : public Function
	{
	public:
		effortFunctionWayPoint(PTPWayPointOptimization* PTPOptimizer) : _PTPOptimizer(PTPOptimizer) {}

		VectorX func(const VectorX& params) const;
		MatrixX Jacobian(const VectorX& params) const;

		PTPWayPointOptimization* _PTPOptimizer;
	};

	class NonlinearInequalityConstraintWayPoint : public Function
	{
	public:
		NonlinearInequalityConstraintWayPoint(PTPWayPointOptimization* PTPOptimizer) : _PTPOptimizer(PTPOptimizer) {}

		VectorX func(const VectorX& params) const;
		MatrixX Jacobian(const VectorX& params) const;

		PTPWayPointOptimization* _PTPOptimizer;
	};

	//class energyLossFunction : public Function
	//{
	//	energyLossFunction() {}

	//	VectorX func(const VectorX& x) const;
	//	MatrixX Jacobian(const VectorX& x) const;
	//};

	//class equalityConstraint : public Function
	//{
	//	equalityConstraint() {}

	//	VectorX func(const VectorX& x) const;
	//	MatrixX Jacobian(const VectorX& x) const;
	//};
}