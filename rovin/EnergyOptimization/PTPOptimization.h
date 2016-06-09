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
#include <rovin/Math/GCMMAOptimization.h>
#include <vector>


namespace rovin {

	class PTPOptimization;
	class sharedResource;

	class effortFunction;
	class energyLossFunction;
	class InequalityConstraint;
	
	typedef std::shared_ptr<sharedResource> sharedResourcePtr;

	enum ObjectiveFunctionType { effort, energyloss };
	enum OptimizationType { nlopt, GCMMA, GCMMA_TR, GCMMA_GD };

	// 제로 속도에서 제로 속도로 가는거..?!!!

	class PTPOptimization
	{
		friend class sharedResource;
		friend class effortFunction;
		friend class NonlinearInequalityConstraint;

	public:
		

	public:
		NonlinearOptimization _optimizer;
		GCMMAOptimization * _GCMMAoptimizer;
//		GCMMA_PDIPM _GCMMAoptimizer;

		OptimizationType _optType;
		ObjectiveFunctionType _objectiveType;

		sharedResourcePtr _shared;
		FunctionPtr _objectFunc;
		//equalityConstraint _equlfunc;

		FunctionPtr _IneqFunc;
		FunctionPtr _linearIneqFunc;
		FunctionPtr _nonlinearIneqFunc;

		SerialOpenChainPtr _soc; ///< Serial open chain robot
		
		std::vector<bool> _optJoint;
		unsigned int _numOfOptJoint;
		std::vector<unsigned int> _optJointIdx;
		std::vector<unsigned int> _noptJointIdx;

		unsigned int _orderOfBSpline; ///< order of B spline
		unsigned int _numOfOptCP; ///< number of B spline control points for optimization
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
		
		// initial value
		VectorX initX;

	public:
		PTPOptimization(const SerialOpenChainPtr& soc, const std::vector<bool>& optJoint, const unsigned int orderOfBSpline,
			const unsigned int numOfOptCP, const unsigned int numOfGQSample, const Real tf, const StatePtr& initialState, const StatePtr& finalState);
		PTPOptimization(const SerialOpenChainPtr& soc, const std::vector<bool>& optJoint, const unsigned int orderOfBSpline,
			const unsigned int numOfOptCP, const unsigned int numOfGQSample, const Real tf, const StatePtr& initialState, const StatePtr& finalState, OptimizationType optType);
		PTPOptimization(const SerialOpenChainPtr& soc, const std::vector<bool>& optJoint, const unsigned int orderOfBSpline,
			const unsigned int numOfOptCP, const unsigned int numOfGQSample, const Real tf, const StatePtr& initialState, const StatePtr& finalState, OptimizationType optType, ObjectiveFunctionType objectiveType);
		~PTPOptimization() { delete _GCMMAoptimizer; }

		void makeNonOptJointCP();

		// B spline function
		void makeBSplineKnot();
		void makeBoundaryCondition();

		// Optimization function
		void makeObjectiveFunction();
		void makeIneqConstraintFunction_nlopt();
		void makeIneqConstraintFunction_MMA();
		void generateTrajectory();
	};

	// calculate inverse 
	class sharedResource
	{
		friend class PTPOptimization;
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
		sharedResource(PTPOptimization* PTPOptimizer);

		void makeBSpline(const VectorX& params);

		void update(const VectorX& params);

		//
		const std::vector<VectorX>& gettau(const VectorX& params);
		const std::vector<MatrixX>& getdtaudp(const VectorX& params);

		const std::vector<StatePtr> getstate() const { return _state; }
		const std::vector<MatrixX> getdqddotdp() const { return _dqddotdp; }
		const std::vector<MatrixX> getdqdotdp() const { return _dqdotdp; }

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
	};



	class effortFunction : public Function
	{
	public:
		effortFunction(PTPOptimization* PTPOptimizer) : _PTPOptimizer(PTPOptimizer) {}

		VectorX func(const VectorX& params) const;
		MatrixX Jacobian(const VectorX& params) const;

#ifdef STRATEGY_SCALE
		VectorX func(const VectorX& params, const Real& scObjFunc, const VectorX& scParams) 
		{
			VectorX tmpVec = params;
			for (int j = 0; j < tmpVec.size(); j++)
				tmpVec(j) /= scParams(j);
			VectorX fval = func(tmpVec);
			fval(0) *= scObjFunc;
			return fval;
		}
		MatrixX Jacobian(const VectorX& params, const Real& scObjFunc, const VectorX& scParams) 
		{
			VectorX tmpVec = params;
			for (int j = 0; j < tmpVec.size(); j++)
				tmpVec(j) /= scParams(j);
			MatrixX jacobian = Jacobian(tmpVec);
			for (int j = 0; j < tmpVec.size(); j++)
				jacobian(0, j) *= scObjFunc / scParams(j);
			return jacobian;
		}
#endif

		PTPOptimization* _PTPOptimizer;
	};

	class energyLossFunction : public Function
	{
	public:
		energyLossFunction(PTPOptimization* PTPOptimizer) : _PTPOptimizer(PTPOptimizer) {}

		VectorX func(const VectorX& params) const;
		MatrixX Jacobian(const VectorX& params) const;

#ifdef STRATEGY_SCALE
		VectorX func(const VectorX& params, const Real& scObjFunc, const VectorX& scParams) 
		{
			VectorX tmpVec = params;
			for (int j = 0; j < tmpVec.size(); j++)
				tmpVec(j) /= scParams(j);
			VectorX fval = func(tmpVec);
			fval(0) *= scObjFunc;
			return fval;
		}
		MatrixX Jacobian(const VectorX& params, const Real& scObjFunc, const VectorX& scParams) 
		{
			VectorX tmpVec = params;
			for (int j = 0; j < tmpVec.size(); j++)
				tmpVec(j) /= scParams(j);
			MatrixX jacobian = Jacobian(tmpVec);
			for (int j = 0; j < tmpVec.size(); j++)
				jacobian(0, j) *= scObjFunc / scParams(j);
			return jacobian;
		}
#endif

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

#ifdef STRATEGY_SCALE

	class scaleAugmentedFunction : public Function
	{
	public:
		scaleAugmentedFunction() {}

		VectorX func(const VectorX& x) const
		{
			std::vector<VectorX> fvalList(_functionList.size());
			unsigned int dimension = 0;
			for (unsigned int i = 0; i < _functionList.size(); i++)
			{
				fvalList[i] = _functionList[i]->func(x);
				dimension += fvalList[i].size();
			}
			VectorX f(dimension);
			for (unsigned int i = 0, idx = 0; i < _functionList.size(); i++)
			{
				f.block(idx, 0, fvalList[i].size(), 1) = fvalList[i];
				idx += fvalList[i].size();
			}
			return f;
		}
		MatrixX Jacobian(const VectorX& x) const
		{
			std::vector<MatrixX> jacobianList(_functionList.size());
			unsigned int dimension = 0;
			for (unsigned int i = 0; i < _functionList.size(); i++)
			{
				jacobianList[i] = _functionList[i]->Jacobian(x);
				dimension += jacobianList[i].rows();
			}
			MatrixX jacobian(dimension, x.size());
			for (unsigned int i = 0, idx = 0; i < _functionList.size(); i++)
			{
				jacobian.block(idx, 0, jacobianList[i].rows(), x.size()) = jacobianList[i];
				idx += jacobianList[i].rows();
			}
			return jacobian;
		}

		void addFunction(const FunctionPtr& function)
		{
			_functionList.push_back(function);
		}


		VectorX func(const VectorX& params, const VectorX& scIneqFunc, const VectorX& scParams) 
		{
			VectorX tmpVec = params;
			for (int j = 0; j < tmpVec.rows(); j++)
				tmpVec(j) /= scParams(j);
			VectorX fval = func(tmpVec);
			for (int i = 0; i < fval.size(); i++)
				fval(i) *= scIneqFunc(i);
			return fval;
		}
		MatrixX Jacobian(const VectorX& params, const VectorX& scIneqFunc, const VectorX& scParams) 
		{
			VectorX tmpVec = params;
			for (int j = 0; j < tmpVec.rows(); j++)
				tmpVec(j) /= scParams(j);
			MatrixX jacobian = Jacobian(tmpVec);
			for (int i = 0; i < jacobian.rows(); i++)
				for (int j = 0; j < jacobian.cols(); j++)
					jacobian(i, j) *= scIneqFunc(i) / scParams(j);
			return jacobian;
		}

	private:
		std::vector<FunctionPtr> _functionList;

	};
#endif


}
