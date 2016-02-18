/*!
 *	\file	PTPOptimization.h
 *	\date	2016.01.22
 *	\author	Youngsuk (crazyhys@gmail.com)
 *	\brief	PTPOptimization class
 *          this class make a energy efficiency optimization path of PToP problem
*/

#pragma once

#include <rovin/Optimizer/NonlinearOptimization.h>

#include <rovin/Dynamics/SerialOpenChain.h>
#include <rovin/Math/GaussianQuadrature.h>
#include <rovin/Math/Interpolation.h>

namespace rovin {

	class PTPOptimization;
	class sharedResource;


	class effortFunction;
	class NonlinearInequalityConstraint;

	class energyLossFunction;
	class EqualityConstraint;
	class InequalityConstraint;
	

	typedef std::shared_ptr<sharedResource> sharedResourcePtr;

	// 제로 속도에서 제로 속도로 가는거..?!!!

	class PTPOptimization
	{
		friend class sharedResource;
		friend class effortFunction;
		friend class NonlinearInequalityConstraint;

	public:
		enum ObjectiveFunctionType { effort, energyloss };
		
	public:
		NonlinearOptimization _optimizer;
		sharedResourcePtr _shared;

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

		unsigned int _orderOfBSpline; ///< order of B spline
		unsigned int _numOfOptCP; ///< number of B spline control points for optimization except contraint control points
		unsigned int _numOfGQSample; ///< number of Gaussian quadrature sampling

		Real _tf; ///< traveling time

		StatePtr _initialState; ///< initial joint position
		StatePtr _finalState; ///< final joint position

		// B spline variables
		VectorX _knot;
		std::vector<VectorX> _initialCP;
		std::vector<VectorX> _finalCP;
		MatrixX _noptJointCP;

		// Gaussian quadrature
		GaussianQuadrature GQ;

	public:
		PTPOptimization(const SerialOpenChainPtr& soc, const std::vector<bool>& optJoint, const unsigned int orderOfBSpline,
			const unsigned int numOfCP, const unsigned int numOfGQSample, const Real tf, const StatePtr& initialState, const StatePtr& finalState);

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
	class sharedResource
	{

	private:

		PTPOptimization* _PTPOptimizer;

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
		sharedResource(PTPOptimization* PTPOptimizer);

		void makeBSpline(const VectorX& params);

		void update(const VectorX& params);

		const std::vector<VectorX>& gettau(const VectorX& params);
		const std::vector<MatrixX>& getdtaudp(const VectorX& params);

	};

	class effortFunction : public Function
	{
	public:
		effortFunction(PTPOptimization* PTPOptimizer) : _PTPOptimizer(PTPOptimizer) {}

		VectorX func(const VectorX& params) const;
		MatrixX Jacobian(const VectorX& params) const;

		PTPOptimization* _PTPOptimizer;
	};

	class NonlinearInequalityConstraint : public Function
	{
	public:
		NonlinearInequalityConstraint(PTPOptimization* PTPOptimizer) : _PTPOptimizer(PTPOptimizer) {}

		VectorX func(const VectorX& params) const;
		MatrixX Jacobian(const VectorX& params) const;

		PTPOptimization* _PTPOptimizer;
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
