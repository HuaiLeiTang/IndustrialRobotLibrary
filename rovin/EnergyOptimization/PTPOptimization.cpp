#include "PTPOptimization.h"

#include <time.h>
#include <memory>

using namespace std;

namespace rovin{

	void PTPOptimization::contructorSetting()
	{
		_dt = (_tf) / (_numOfData - 1);
		Real t = 0.0;
		_tspan.resize(_numOfData);
		for (unsigned int i = 0; i < _numOfData; i++)
		{
			_tspan(i) = t;
			t += _dt;
		}
	}

	PTPOptimization::PTPOptimization(const SerialOpenChainPtr & soc, const vector<bool>& optJoint,
		const unsigned int orderOfBSpline, const unsigned int numOfOptCP, const unsigned int numOfGQSample,
		const Real tf, const StatePtr& initialState, const StatePtr& finalState)
	{
		_soc = soc;
		_optJoint = optJoint;
		_orderOfBSpline = orderOfBSpline;
		_numOfOptCP = numOfOptCP;
		_numOfGQSample = numOfGQSample;
		_tf = tf;
		_initialState = initialState;
		_finalState = finalState;
		_optType = OptimizationType::nlopt;
		_objectiveType = ObjectiveFunctionType::effort;

		_numOfOptJoint = 0;
		for (unsigned int i = 0; i < _soc->getNumOfJoint(); i++)
		{
			if (_optJoint[i])
			{
				_optJointIdx.push_back(i);
				_numOfOptJoint++;
			}
			else
			{
				_noptJointIdx.push_back(i);
			}
		}
		_GCMMAoptimizer = NULL;
		_initialswi = false;
		contructorSetting();
	}


	PTPOptimization::PTPOptimization(const SerialOpenChainPtr& soc, const std::vector<bool>& optJoint, const unsigned int orderOfBSpline,
		const unsigned int numOfOptCP, const unsigned int numOfGQSample, const Real tf, const StatePtr& initialState, const StatePtr& finalState, OptimizationType optType)
	{
		_soc = soc;
		_optJoint = optJoint;
		_orderOfBSpline = orderOfBSpline;
		_numOfOptCP = numOfOptCP;
		_numOfGQSample = numOfGQSample;
		_tf = tf;
		_initialState = initialState;
		_finalState = finalState;

		_numOfOptJoint = 0;
		for (unsigned int i = 0; i < _soc->getNumOfJoint(); i++)
		{
			if (_optJoint[i])
			{
				_optJointIdx.push_back(i);
				_numOfOptJoint++;
			}
			else
			{
				_noptJointIdx.push_back(i);
			}
		}
		_optType = optType;
		_GCMMAoptimizer = NULL;
		_initialswi = false;
		contructorSetting();
	}

	PTPOptimization::PTPOptimization(const SerialOpenChainPtr& soc, const std::vector<bool>& optJoint, const unsigned int orderOfBSpline,
		const unsigned int numOfOptCP, const unsigned int numOfGQSample, const Real tf, const StatePtr& initialState, const StatePtr& finalState, OptimizationType optType, ObjectiveFunctionType objectiveType)
	{
		_soc = soc;
		_optJoint = optJoint;
		_orderOfBSpline = orderOfBSpline;
		_numOfOptCP = numOfOptCP;
		_numOfGQSample = numOfGQSample;
		_tf = tf;
		_initialState = initialState;
		_finalState = finalState;

		_numOfOptJoint = 0;
		for (unsigned int i = 0; i < _soc->getNumOfJoint(); i++)
		{
			if (_optJoint[i])
			{
				_optJointIdx.push_back(i);
				_numOfOptJoint++;
			}
			else
			{
				_noptJointIdx.push_back(i);
			}
		}
		_optType = optType;
		_objectiveType = objectiveType;
		_GCMMAoptimizer = NULL;
		_initialswi = false;
		contructorSetting();
	}

	PTPOptimization::PTPOptimization(const SerialOpenChainPtr& soc, const std::vector<bool>& optJoint, const BSpline<-1, -1, -1>& initialBSpline, const unsigned int numOfGQSample,
		const StatePtr& initialState, const StatePtr& finalState, OptimizationType optType, ObjectiveFunctionType objectiveType)
		: _soc(soc), _optJoint(optJoint), _initialBSpline(initialBSpline), _numOfGQSample(numOfGQSample), _initialState(initialState), 
		_finalState(finalState), _optType(optType), _objectiveType(objectiveType)
	{
		_numOfOptJoint = 0;
		for (unsigned int i = 0; i < _soc->getNumOfJoint(); i++)
		{
			if (_optJoint[i])
			{
				_optJointIdx.push_back(i);
				_numOfOptJoint++;
			}
			else
			{
				_noptJointIdx.push_back(i);
			}
		}
		_GCMMAoptimizer = NULL;

		_orderOfBSpline = _initialBSpline.getOrder();
		_numOfOptCP = _initialBSpline.getControlPoints().cols() - 6;
		_tf = _initialBSpline.getKnots()[_initialBSpline.getKnots().size() - 1];
		_initialswi = true;
		contructorSetting();
	}

	void PTPOptimization::makeBSplineKnot()
	{
		LOGIF(_orderOfBSpline + _numOfOptCP + 6 > 2 * _orderOfBSpline, "The number of control points is not enough.");

		_knot.resize(_orderOfBSpline + _numOfOptCP + 6);
		for (unsigned int i = 0; i < _orderOfBSpline; i++)
		{
			_knot[i] = 0;
			_knot[_knot.size() - i - 1] = _tf;
		}

		Real delta = _tf / (_numOfOptCP + 6 - _orderOfBSpline + 1);
		for (unsigned int i = 0; i < _numOfOptCP + 6 - _orderOfBSpline; i++)
		{
			_knot[_orderOfBSpline + i] = delta * (i + 1);
		}
	}

	void PTPOptimization::makeBoundaryCondition()
	{
		Real delta = _tf / (_numOfOptCP + 6 - _orderOfBSpline + 1);

		_initialCP.resize(3);
		_initialCP[0] = _initialState->getJointStatePos();
		_initialCP[1] = delta / (_orderOfBSpline - 1)*_initialState->getJointStateVel() + _initialCP[0];
		_initialCP[2] = 2 * delta*(delta / (_orderOfBSpline - 1) / (_orderOfBSpline - 2)*_initialState->getJointStateAcc() + _initialCP[1] * (1 / (2 * delta) + 1 / delta) - _initialCP[0] / delta);

		_finalCP.resize(3);
		_finalCP[0] = _finalState->getJointStatePos();
		_finalCP[1] = -delta / (_orderOfBSpline - 1)*_finalState->getJointStateVel() + _finalCP[0];
		_finalCP[2] = 2 * delta*(delta / (_orderOfBSpline - 1) / (_orderOfBSpline - 2)*_finalState->getJointStateAcc() + _finalCP[1] * (1 / (2 * delta) + 1 / delta) - _finalCP[0] / delta);
	}

	void PTPOptimization::makeObjectiveFunction()
	{
		if (_objectiveType == ObjectiveFunctionType::effort)
		{
			_objectFunc = FunctionPtr(new effortFunction(this));
		}
		else if (_objectiveType == ObjectiveFunctionType::energyloss)
		{
			_objectFunc = FunctionPtr(new energyLossFunction(this));
		}
		else if (_objectiveType == ObjectiveFunctionType::acceleration)
		{
			_objectFunc = FunctionPtr(new accelerationFunction(this));
		}
	}

	void PTPOptimization::makeIneqConstraintFunction_nlopt()
	{
		_nonlinearIneqFunc = FunctionPtr(new NonlinearInequalityConstraint(this));

		unsigned int np = _shared->_dPdP.rows(); // _qSpline.getControlPoints().cols() - 6,
		unsigned int nq = _shared->_dQdP.rows(); // _qdotSpline.getControlPoints().cols() - 4
		unsigned int nr = _shared->_dRdP.rows(); // _qddotSpline.getControlPoints().cols() - 2

		MatrixX A((np + nq + nr) * 2 * _numOfOptJoint, _numOfOptCP * _numOfOptJoint);
		VectorX b((np + nq + nr) * 2 * _numOfOptJoint);
		A.setZero();
		b.setZero();
		for (unsigned int i = 0; i < _numOfOptJoint; i++)
		{
			A.block((np + nq + nr) * 2 * i, _numOfOptCP * i, np, _numOfOptCP) = _shared->_dPdP;
			A.block((np + nq + nr) * 2 * i + np, _numOfOptCP * i, np, _numOfOptCP) = -_shared->_dPdP;

			A.block((np + nq + nr) * 2 * i + np * 2, _numOfOptCP * i, nq, _numOfOptCP) = _shared->_dQdP;
			A.block((np + nq + nr) * 2 * i + np * 2 + nq, _numOfOptCP * i, nq, _numOfOptCP) = -_shared->_dQdP;

			A.block((np + nq + nr) * 2 * i + (np + nq) * 2, _numOfOptCP * i, nr, _numOfOptCP) = _shared->_dRdP;
			A.block((np + nq + nr) * 2 * i + (np + nq) * 2 + nr, _numOfOptCP * i, nr, _numOfOptCP) = -_shared->_dRdP;

			b.block((np + nq + nr) * 2 * i, 0, np, 1) = _shared->_P[i] - VectorX::Ones(np)*_soc->getMotorJointPtr(_optJointIdx[i])->getLimitPosUpper();
			b.block((np + nq + nr) * 2 * i + np, 0, np, 1) = -_shared->_P[i] + VectorX::Ones(np)*_soc->getMotorJointPtr(_optJointIdx[i])->getLimitPosLower();

			b.block((np + nq + nr) * 2 * i + np * 2, 0, nq, 1) = _shared->_Q[i] -VectorX::Ones(nq)*_soc->getMotorJointPtr(_optJointIdx[i])->getLimitVelUpper();
			b.block((np + nq + nr) * 2 * i + np * 2 + nq, 0, nq, 1) = -_shared->_Q[i] + VectorX::Ones(nq)*_soc->getMotorJointPtr(_optJointIdx[i])->getLimitVelLower();

			b.block((np + nq + nr) * 2 * i + (np + nq) * 2, 0, nr, 1) = _shared->_R[i] -VectorX::Ones(nr)*_soc->getMotorJointPtr(_optJointIdx[i])->getLimitAccUpper();
			b.block((np + nq + nr) * 2 * i + (np + nq) * 2 + nr, 0, nr, 1) = -_shared->_R[i] + VectorX::Ones(nr)*_soc->getMotorJointPtr(_optJointIdx[i])->getLimitAccLower();
		}
		_linearIneqFunc = FunctionPtr(new AffineFunction(A, b));

		_IneqFunc = FunctionPtr(new AugmentedFunction());
		static_pointer_cast<AugmentedFunction>(_IneqFunc)->addFunction(_nonlinearIneqFunc);
		static_pointer_cast<AugmentedFunction>(_IneqFunc)->addFunction(_linearIneqFunc);

	}

	void PTPOptimization::makeIneqConstraintFunction_MMA()
	{
		_nonlinearIneqFunc = FunctionPtr(new NonlinearInequalityConstraint(this));

		unsigned int nq = _shared->_dQdP.rows(); // _qdotSpline.getControlPoints().cols() - 4
		unsigned int nr = _shared->_dRdP.rows(); // _qddotSpline.getControlPoints().cols() - 2

		MatrixX A((nq + nr) * 2 * _numOfOptJoint, _numOfOptCP * _numOfOptJoint);
		VectorX b((nq + nr) * 2 * _numOfOptJoint);
		A.setZero();
		b.setZero();
		for (unsigned int i = 0; i < _numOfOptJoint; i++)
		{
			A.block((nq + nr) * 2 * i, _numOfOptCP * i, nq, _numOfOptCP) = _shared->_dQdP;
			A.block((nq + nr) * 2 * i + nq, _numOfOptCP * i, nq, _numOfOptCP) = -_shared->_dQdP;

			A.block((nq + nr) * 2 * i + (nq)* 2, _numOfOptCP * i, nr, _numOfOptCP) = _shared->_dRdP;
			A.block((nq + nr) * 2 * i + (nq)* 2 + nr, _numOfOptCP * i, nr, _numOfOptCP) = -_shared->_dRdP;

			b.block((nq + nr) * 2 * i, 0, nq, 1) = _shared->_Q[i] - VectorX::Ones(nq)*_soc->getMotorJointPtr(_optJointIdx[i])->getLimitVelUpper();
			b.block((nq + nr) * 2 * i + nq, 0, nq, 1) = -_shared->_Q[i] + VectorX::Ones(nq)*_soc->getMotorJointPtr(_optJointIdx[i])->getLimitVelLower();

			b.block((nq + nr) * 2 * i + (nq)* 2, 0, nr, 1) = _shared->_R[i] - VectorX::Ones(nr)*_soc->getMotorJointPtr(_optJointIdx[i])->getLimitAccUpper();
			b.block((nq + nr) * 2 * i + (nq)* 2 + nr, 0, nr, 1) = -_shared->_R[i] + VectorX::Ones(nr)*_soc->getMotorJointPtr(_optJointIdx[i])->getLimitAccLower();
		}
		_linearIneqFunc = FunctionPtr(new AffineFunction(A, b));

		_IneqFunc = FunctionPtr(new AugmentedFunction());
		static_pointer_cast<AugmentedFunction>(_IneqFunc)->addFunction(_nonlinearIneqFunc);
		static_pointer_cast<AugmentedFunction>(_IneqFunc)->addFunction(_linearIneqFunc);
	}

	void PTPOptimization::makeNonOptJointCP()
	{
		_noptJointCP.resize(_soc->getNumOfJoint() - _numOfOptJoint, _numOfOptCP);
		for (unsigned int i = 0; i < _soc->getNumOfJoint() - _numOfOptJoint; i++)
		{
			for (unsigned int j = 0; j < _numOfOptCP; j++)
			{
				_noptJointCP(i, j) = (_finalCP[2](_noptJointIdx[i]) - _initialCP[2](_noptJointIdx[i])) / (_numOfOptCP + 1) * (j + 1) + _initialCP[2](_noptJointIdx[i]);
			}
		}
	}

	void PTPOptimization::generateTrajectory()
	{
		if(!_initialswi)
			makeBSplineKnot();
		else
		{
			_knot = _initialBSpline.getKnots();
		}
		

		LOG("Complete making B-Spline knots.");
		//cout << "knot : " << _knot.transpose() << endl;
		if (!_initialswi)
			makeBoundaryCondition();
		else
		{
			_initialCP.resize(3); _finalCP.resize(3);
			for (unsigned int i = 0; i < 3; i++)
			{
				_initialCP[i] = _initialBSpline.getControlPoints().col(i);
				_finalCP[i] = _initialBSpline.getControlPoints().col(_initialBSpline.getControlPoints().cols() - 1 - i);
			}
		}

		LOG("Complete making boundary conditions.");
		//cout << "Boundary Conditions : " << endl;
		//for (unsigned int i = 0; i < 3; i++) cout << _initialCP[i].transpose() << endl;
		//for (unsigned int i = 0; i < 3; i++) cout << _finalCP[i].transpose() << endl; 

		GQ = GaussianQuadrature(_numOfGQSample, 0.0, _tf);
		LOG("Gaussian quadrature ready.");
		//cout << "Gaussian Qudrature ti: " << GQ.getQueryPoints().transpose() << endl;
		//cout << "Gaussian Qudrature wi: " << GQ.getWeights().transpose() << endl;

		_shared = sharedResourcePtr(new sharedResource(this));
		LOG("Complete shared initialization.");
		//cout << _shared->_dPdP << endl;
		//cout << _shared->_dQdP << endl;
		//cout << _shared->_dRdP << endl;

		if (!_initialswi)
			makeNonOptJointCP();
		else
			_noptJointCP = _initialBSpline.getControlPoints().block(3, 3, 3, _numOfOptCP);

		makeObjectiveFunction();
		if (_optType == OptimizationType::nlopt)
		{
			makeIneqConstraintFunction_nlopt();
		}
		else if ( (_optType == OptimizationType::GCMMA) || (_optType == OptimizationType::GCMMA_TR) || (_optType == OptimizationType::GCMMA_GD))
		{
			makeIneqConstraintFunction_nlopt();
			//makeIneqConstraintFunction_MMA();
		}
		
		LOG("Optimization ready.");
		//cout << "Nopt Joint CP: " << endl;
		//cout << _noptJointCP << endl;
		//cout << "LinearInequalityCondition(A): " << endl << static_pointer_cast<AffineFunction>(_linearIneqFunc)->getA() << endl;
		//cout << "LinearInequalityCondition(b): " << endl << static_pointer_cast<AffineFunction>(_linearIneqFunc)->getb() << endl;

		// Initial Guess
		initX.resize(_numOfOptJoint * _numOfOptCP);

		if (!_initialswi)
		{
			srand((unsigned)time(NULL));
			for (unsigned int i = 0; i < _numOfOptJoint; i++)
			{
				for (unsigned int j = 0; j < _numOfOptCP; j++)
				{
					initX(_numOfOptCP * i + j) = (_finalCP[2](_optJointIdx[i]) - _initialCP[2](_optJointIdx[i])) / (_numOfOptCP + 1) * (j + 1) + _initialCP[2](_optJointIdx[i]);
					//initX(_numOfOptCP * i + j) = _soc->getMotorJointPtr(i)->getLimitPosLower() + (_soc->getMotorJointPtr(i)->getLimitPosUpper() - _soc->getMotorJointPtr(i)->getLimitPosLower())*((double)rand() / 32767.0);
				}
			}
		}
		else
		{
			for (unsigned int i = 0; i < _numOfOptJoint; i++)
			{
				for (unsigned int j = 0; j < _numOfOptCP; j++)
				{
					initX(_numOfOptCP * i + j) = _initialBSpline.getControlPoints()(i, j + 3);
				}
			}
		}

		FunctionPtr f = FunctionPtr(new energyLossFunction(this));

		//cout << "_numOfOptJoint : " << _numOfOptJoint << endl;
		//cout << "_orderOfBSpline : " << _orderOfBSpline << endl;
		//cout << "_numOfOptCP : " << _numOfOptCP << endl;
		//cout << "_knots : " << _knot << endl;
		//cout << "initX" << endl << initX << endl;
		//cout << "initialCP" << endl;
		//for (unsigned int i = 0; i < _initialCP.size(); i++)
		//	cout << _initialCP[i] << endl << endl;
		//cout << endl;
		//cout << "_finalCP" << endl;
		//for (unsigned int i = 0; i < _finalCP.size(); i++)
		//	cout << _finalCP[i] << endl << endl;
		//cout << endl;
		//cout << "_noptJointCP" << endl << _noptJointCP << endl;
		//cout << "tf : " << _tf << endl;
		//cout << "_optJointIdx" << endl;
		//for (int i = 0; i < _optJointIdx.size(); i++)
		//	cout << _optJointIdx[i] << '\t';
		//cout << endl;
		//cout << "_noptJointIdx" << endl;
		//for (int i = 0; i < _noptJointIdx.size(); i++)
		//	cout << _noptJointIdx[i] << '\t';
		//cout << endl;

		///////////
		//unsigned int dataNum = initM.cols();
		//Real t = 0.0;
		//Real dt = (_tf - t) / (dataNum - 1);
		//VectorX tt(dataNum);
		//for (unsigned int i = 0; i < dataNum; i++)
		//{
		//	tt(i) = t;
		//	t += dt;
		//}

		//BSpline<-1, -1, -1> q = BSplineFitting(initM, _orderOfBSpline, _numOfOptCP + 6, tt);
		//MatrixX cptmp = q.getControlPoints();
		//initX << cptmp(0, 3), cptmp(0, 4), cptmp(0, 5), cptmp(0, 6), cptmp(1, 3), cptmp(1, 4), cptmp(1, 5), cptmp(1, 6), cptmp(2, 3), cptmp(2, 4), cptmp(2, 5), cptmp(2, 6);

		//cout << "initCP" << endl;
		//for (int i = 0; i < _initialCP.size(); i++)
		//	cout << _initialCP[i] << endl;
		//cout << "finalCP" << endl;
		//for (int i = 0; i < _finalCP.size(); i++)
		//	cout << _finalCP[i] << endl;
		//cout << "_noptJointCP" << endl;
		//cout << _noptJointCP << endl;

		//initX << 0.556274, 0.11917, -0.328995, -0.814881, 0.603234, 0.0932657, -0.130955, -0.064882, 1.14852, 0.774952, 0.2058, -0.233091;
		//initX << 0.419949, 0.275863, -0.263421, -0.847695, -0.218918, -0.226281, -0.224308, -0.177036, 1.22172, 1.20281, 0.169464, -0.362267;
		LOG("Initial guess ready.");
		//cout << "Initial X: " << endl;%
		//cout << initX << endl;
		cout << "f(X): " << endl;
		cout << _objectFunc->func(initX) << endl;
		cout << "energy loss before optimization : " << f->func(initX) << endl;
		//cout << "Inequality(X): " << endl;
		//cout << _IneqFunc->func(initX) << endl;

		if (_optType == OptimizationType::nlopt)
		{
			_optimizer.setObjectiveFunction(_objectFunc);
			_optimizer.setInequalityConstraint(_IneqFunc);
			LOG("Start optimization.");
			clock_t time = clock();
			//for (int i = 0; i < 100; i++)
				_optimizer.solve(initX);
			LOG("Finish optimization.");
			cout << "------------------------------------" << endl;
			cout << "computation time : " << (clock() - time) << endl << endl;
			cout << "X : " << endl << _optimizer.resultX << endl << endl;
			//_shared->makeBSpline(_optimizer.resultX);
			//cout << _knot << endl;
			//cout << "control points" << endl << _shared->_qSpline.getControlPoints() << endl << endl;
			//cout << "Inequality : " << _IneqFunc->func(_optimizer.resultX) << endl;
			cout << "Value of objective function : " << _optimizer.resultFunc << endl << endl;
			cout << "energy loss after optimization : " << f->func(_optimizer.resultX) << endl;
		}
		else if ((_optType == OptimizationType::GCMMA) || (_optType == OptimizationType::GCMMA_TR) || (_optType == OptimizationType::GCMMA_GD))
		{
			VectorX minX(initX.size()), maxX(initX.size());
			for (unsigned int iii = 0; iii < _numOfOptJoint; iii++)
			{
				for (unsigned int jjj = 0; jjj < _numOfOptCP; jjj++)
				{
					minX(iii * _numOfOptCP + jjj) = _soc->getMotorJointPtr(_optJointIdx[iii])->getLimitPosLower();
					maxX(iii * _numOfOptCP + jjj) = _soc->getMotorJointPtr(_optJointIdx[iii])->getLimitPosUpper();
				}
			}

			if (_optType == OptimizationType::GCMMA)
				_GCMMAoptimizer = new GCMMA_PDIPM(initX.size(), _IneqFunc->func(initX).size());
			else if (_optType == OptimizationType::GCMMA_TR)
				_GCMMAoptimizer = new GCMMA_TRM(initX.size(), _IneqFunc->func(initX).size());
			else if (_optType == OptimizationType::GCMMA_GD)
				_GCMMAoptimizer = new GCMMA_GDM(initX.size(), _IneqFunc->func(initX).size());

			_GCMMAoptimizer->setObjectiveFunction(_objectFunc);
			_GCMMAoptimizer->setInequalityConstraint(_IneqFunc);
			_GCMMAoptimizer->setMinMax(minX, maxX);

			//VectorX fvali = _objectFunc->func(initX);
			//VectorX InequalVal = _IneqFunc->func(initX);

			LOG("Start optimization.");
			clock_t time = clock();
			GCMMAReturnFlag retFlag;
			retFlag = _GCMMAoptimizer->solve(initX);
			LOG("Finish optimization.");
			cout << "------------------------------------" << endl;
			cout << "computation time : " << (clock() - time) << endl << endl;

			displayGCMMAResult(retFlag);
			cout << "X : " << endl << _GCMMAoptimizer->getResultX() << endl << endl;
			
			//VectorX X = _GCMMAoptimizer->getResultX();
			//cout << "func(X) : " << _GCMMAoptimizer->_objectFunc->func(X) << endl;
			//cout << "ineqCons : " << _GCMMAoptimizer->_ineqConstraint->func(X) << endl;

			cout << "control points" << endl << _shared->_qSpline.getControlPoints() << endl << endl;
			cout << "Value of objective function : " << _GCMMAoptimizer->getResultFunc() << endl << endl;

			//if (_optType == OptimizationType::GCMMA)
			//	cout << "_suby : " << endl << _GCMMAoptimizer._suby << endl << endl;
			//else if (_optType == OptimizationType::GCMMA_TR)
			//	cout << "_suby : " << endl << _GCMMAoptimizer._TR_suby << endl << endl;
		}
	}
	
	// class sharedResource

	sharedResource::sharedResource(PTPOptimization* PTPOptimizer)
	{
		//////////////////////////////// choi /////////////////////////////////////////
		//_PTPOptimizer = PTPOptimizer;
		//_state.resize(_PTPOptimizer->_numOfGQSample);
		//for (unsigned int i = 0; i < _state.size(); i++)
		//{
		//	_state[i] = _PTPOptimizer->_soc->makeState();
		//}

		//// Calculate _dqdp,dqdotdp, dqddotdq
		//// Calculate _dPdP,_dQdP, _dRdP
		//MatrixX cp(1, _PTPOptimizer->_numOfOptCP + 6);
		//bool checkMatrixSize = false;
		//_dqdp.resize(_PTPOptimizer->_numOfGQSample, MatrixX::Zero(_PTPOptimizer->_soc->getNumOfJoint(), _PTPOptimizer->_numOfOptCP * _PTPOptimizer->_numOfOptJoint));
		//_dqdotdp.resize(_PTPOptimizer->_numOfGQSample, MatrixX::Zero(_PTPOptimizer->_soc->getNumOfJoint(), _PTPOptimizer->_numOfOptCP * _PTPOptimizer->_numOfOptJoint));
		//_dqddotdp.resize(_PTPOptimizer->_numOfGQSample, MatrixX::Zero(_PTPOptimizer->_soc->getNumOfJoint(), _PTPOptimizer->_numOfOptCP * _PTPOptimizer->_numOfOptJoint));
		//for (unsigned int i = 0; i < _PTPOptimizer->_numOfOptCP; i++)
		//{
		//	cp.setZero();
		//	cp(0, 3 + i) = 1.0;
		//	_qSpline = BSpline<-1, -1, -1>(_PTPOptimizer->_knot, cp);
		//	_qdotSpline = _qSpline.derivative();
		//	_qddotSpline = _qdotSpline.derivative();

		//	if (!checkMatrixSize)
		//	{
		//		_dPdP.resize(_qSpline.getControlPoints().cols() - 6, _PTPOptimizer->_numOfOptCP);
		//		_dQdP.resize(_qdotSpline.getControlPoints().cols() - 4, _PTPOptimizer->_numOfOptCP);
		//		_dRdP.resize(_qddotSpline.getControlPoints().cols() - 2, _PTPOptimizer->_numOfOptCP);
		//		checkMatrixSize = true;
		//	}
		//	_dPdP.col(i) = _qSpline.getControlPoints().block(0, 3, 1, _dPdP.rows()).transpose();
		//	_dQdP.col(i) = _qdotSpline.getControlPoints().block(0, 2, 1, _dQdP.rows()).transpose();
		//	_dRdP.col(i) = _qddotSpline.getControlPoints().block(0, 1, 1, _dRdP.rows()).transpose();

		//	for (unsigned int j = 0; j < _PTPOptimizer->_numOfGQSample; j++)
		//	{
		//		VectorX dqdp = _qSpline(_PTPOptimizer->GQ.getQueryPoints()[j]);
		//		VectorX dqdotdp = _qdotSpline(_PTPOptimizer->GQ.getQueryPoints()[j]);
		//		VectorX dqddotdp = _qddotSpline(_PTPOptimizer->GQ.getQueryPoints()[j]);

		//		for (unsigned int k = 0; k < _PTPOptimizer->_numOfOptJoint; k++)
		//		{
		//			_dqdp[j](_PTPOptimizer->_optJointIdx[k], _PTPOptimizer->_numOfOptCP*k + i) = dqdp[0];
		//			_dqdotdp[j](_PTPOptimizer->_optJointIdx[k], _PTPOptimizer->_numOfOptCP*k + i) = dqdotdp[0];
		//			_dqddotdp[j](_PTPOptimizer->_optJointIdx[k], _PTPOptimizer->_numOfOptCP*k + i) = dqddotdp[0];
		//		}
		//	}
		//}

		////cout << "_dPdP size : " << _dPdP.rows() << '\t' << _dPdP.cols() << endl;
		////cout << "_dQdP size : " << _dQdP.rows() << '\t' << _dQdP.cols() << endl;
		////cout << "_dRdP size : " << _dRdP.rows() << '\t' << _dRdP.cols() << endl;
		////cout << "_dPdP" << endl << _dPdP << endl;
		////cout << "_dQdP" << endl << _dQdP << endl;
		////cout << "_dRdP" << endl << _dRdP << endl;

		//_P.resize(_PTPOptimizer->_numOfOptJoint);
		//_Q.resize(_PTPOptimizer->_numOfOptJoint);
		//_R.resize(_PTPOptimizer->_numOfOptJoint);
		//for (unsigned int i = 0; i < _PTPOptimizer->_numOfOptJoint; i++)
		//{
		//	cp.setZero();
		//	cp(0) = _PTPOptimizer->_initialCP[0](_PTPOptimizer->_optJointIdx[i]);
		//	cp(1) = _PTPOptimizer->_initialCP[1](_PTPOptimizer->_optJointIdx[i]);
		//	cp(2) = _PTPOptimizer->_initialCP[2](_PTPOptimizer->_optJointIdx[i]);
		//	cp(_PTPOptimizer->_numOfOptCP + 3) = _PTPOptimizer->_finalCP[2](_PTPOptimizer->_optJointIdx[i]);
		//	cp(_PTPOptimizer->_numOfOptCP + 4) = _PTPOptimizer->_finalCP[1](_PTPOptimizer->_optJointIdx[i]);
		//	cp(_PTPOptimizer->_numOfOptCP + 5) = _PTPOptimizer->_finalCP[0](_PTPOptimizer->_optJointIdx[i]);
		//	//cout << cp << endl;
		//	_qSpline = BSpline<-1, -1, -1>(_PTPOptimizer->_knot, cp);
		//	_qdotSpline = _qSpline.derivative();
		//	_qddotSpline = _qdotSpline.derivative();
		//	_P[i] = _qSpline.getControlPoints().block(0, 3, 1, _dPdP.rows()).transpose();
		//	_Q[i] = _qdotSpline.getControlPoints().block(0, 2, 1, _dQdP.rows()).transpose();
		//	_R[i] = _qddotSpline.getControlPoints().block(0, 1, 1, _dRdP.rows()).transpose();

		//	//cout << "_P[" << i << "]" << endl << _P[i] << endl;
		//	//cout << "_Q[" << i << "]" << endl << _Q[i] << endl;
		//	//cout << "_R[" << i << "]" << endl << _R[i] << endl;
		//}
		/////////////////////////////////////////////////////////////////////////////////////

		//////////////////////////////// YS /////////////////////////////////////////
		_PTPOptimizer = PTPOptimizer;
		_state.resize(_PTPOptimizer->_numOfData);
		for (unsigned int i = 0; i < _state.size(); i++)
		{
			_state[i] = _PTPOptimizer->_soc->makeState();
		}

		// Calculate _dqdp,dqdotdp, dqddotdq
		// Calculate _dPdP,_dQdP, _dRdP
		MatrixX cp(1, _PTPOptimizer->_numOfOptCP + 6);
		bool checkMatrixSize = false;
		_dqdp.resize(_PTPOptimizer->_numOfData, MatrixX::Zero(_PTPOptimizer->_soc->getNumOfJoint(), _PTPOptimizer->_numOfOptCP * _PTPOptimizer->_numOfOptJoint));
		_dqdotdp.resize(_PTPOptimizer->_numOfData, MatrixX::Zero(_PTPOptimizer->_soc->getNumOfJoint(), _PTPOptimizer->_numOfOptCP * _PTPOptimizer->_numOfOptJoint));
		_dqddotdp.resize(_PTPOptimizer->_numOfData, MatrixX::Zero(_PTPOptimizer->_soc->getNumOfJoint(), _PTPOptimizer->_numOfOptCP * _PTPOptimizer->_numOfOptJoint));
		for (unsigned int i = 0; i < _PTPOptimizer->_numOfOptCP; i++)
		{
			cp.setZero();
			cp(0, 3 + i) = 1.0;
			_qSpline = BSpline<-1, -1, -1>(_PTPOptimizer->_knot, cp);
			_qdotSpline = _qSpline.derivative();
			_qddotSpline = _qdotSpline.derivative();

			if (!checkMatrixSize)
			{
				_dPdP.resize(_qSpline.getControlPoints().cols() - 6, _PTPOptimizer->_numOfOptCP);
				_dQdP.resize(_qdotSpline.getControlPoints().cols() - 4, _PTPOptimizer->_numOfOptCP);
				_dRdP.resize(_qddotSpline.getControlPoints().cols() - 2, _PTPOptimizer->_numOfOptCP);
				checkMatrixSize = true;
			}
			_dPdP.col(i) = _qSpline.getControlPoints().block(0, 3, 1, _dPdP.rows()).transpose();
			_dQdP.col(i) = _qdotSpline.getControlPoints().block(0, 2, 1, _dQdP.rows()).transpose();
			_dRdP.col(i) = _qddotSpline.getControlPoints().block(0, 1, 1, _dRdP.rows()).transpose();

			for (unsigned int j = 0; j < _PTPOptimizer->_numOfData; j++)
			{
				VectorX dqdp = _qSpline(_PTPOptimizer->GQ.getQueryPoints()[j]);
				VectorX dqdotdp = _qdotSpline(_PTPOptimizer->GQ.getQueryPoints()[j]);
				VectorX dqddotdp = _qddotSpline(_PTPOptimizer->GQ.getQueryPoints()[j]);

				for (unsigned int k = 0; k < _PTPOptimizer->_numOfOptJoint; k++)
				{
					_dqdp[j](_PTPOptimizer->_optJointIdx[k], _PTPOptimizer->_numOfOptCP*k + i) = dqdp[0];
					_dqdotdp[j](_PTPOptimizer->_optJointIdx[k], _PTPOptimizer->_numOfOptCP*k + i) = dqdotdp[0];
					_dqddotdp[j](_PTPOptimizer->_optJointIdx[k], _PTPOptimizer->_numOfOptCP*k + i) = dqddotdp[0];
				}
			}
		}
		_P.resize(_PTPOptimizer->_numOfOptJoint);
		_Q.resize(_PTPOptimizer->_numOfOptJoint);
		_R.resize(_PTPOptimizer->_numOfOptJoint);
		for (unsigned int i = 0; i < _PTPOptimizer->_numOfOptJoint; i++)
		{
			cp.setZero();
			cp(0) = _PTPOptimizer->_initialCP[0](_PTPOptimizer->_optJointIdx[i]);
			cp(1) = _PTPOptimizer->_initialCP[1](_PTPOptimizer->_optJointIdx[i]);
			cp(2) = _PTPOptimizer->_initialCP[2](_PTPOptimizer->_optJointIdx[i]);
			cp(_PTPOptimizer->_numOfOptCP + 3) = _PTPOptimizer->_finalCP[2](_PTPOptimizer->_optJointIdx[i]);
			cp(_PTPOptimizer->_numOfOptCP + 4) = _PTPOptimizer->_finalCP[1](_PTPOptimizer->_optJointIdx[i]);
			cp(_PTPOptimizer->_numOfOptCP + 5) = _PTPOptimizer->_finalCP[0](_PTPOptimizer->_optJointIdx[i]);
			//cout << cp << endl;
			_qSpline = BSpline<-1, -1, -1>(_PTPOptimizer->_knot, cp);
			_qdotSpline = _qSpline.derivative();
			_qddotSpline = _qdotSpline.derivative();
			_P[i] = _qSpline.getControlPoints().block(0, 3, 1, _dPdP.rows()).transpose();
			_Q[i] = _qdotSpline.getControlPoints().block(0, 2, 1, _dQdP.rows()).transpose();
			_R[i] = _qddotSpline.getControlPoints().block(0, 1, 1, _dRdP.rows()).transpose();

			//cout << "_P[" << i << "]" << endl << _P[i] << endl;
			//cout << "_Q[" << i << "]" << endl << _Q[i] << endl;
			//cout << "_R[" << i << "]" << endl << _R[i] << endl;
		}
	}

	void sharedResource::makeBSpline(const VectorX& params)
	{
		MatrixX cp(_PTPOptimizer->_soc->getNumOfJoint(), _PTPOptimizer->_numOfOptCP + 6);

		// Boundary
		cp.col(0) = _PTPOptimizer->_initialCP[0];
		cp.col(1) = _PTPOptimizer->_initialCP[1];
		cp.col(2) = _PTPOptimizer->_initialCP[2];
		cp.col(cp.cols() - 1) = _PTPOptimizer->_finalCP[0];
		cp.col(cp.cols() - 2) = _PTPOptimizer->_finalCP[1];
		cp.col(cp.cols() - 3) = _PTPOptimizer->_finalCP[2];

		// Optimized Joint
		for (unsigned int i = 0; i < _PTPOptimizer->_numOfOptJoint; i++)
		{
			for (unsigned int j = 0; j < _PTPOptimizer->_numOfOptCP; j++)
			{
				cp(_PTPOptimizer->_optJointIdx[i], 3 + j) = params(_PTPOptimizer->_numOfOptCP*i + j);
			}
		}

		// Non-optimized Joint
		for (unsigned int i = 0; i < _PTPOptimizer->_soc->getNumOfJoint() - _PTPOptimizer->_numOfOptJoint; i++)
		{
			for (unsigned int j = 0; j < _PTPOptimizer->_numOfOptCP; j++)
			{
				cp(_PTPOptimizer->_noptJointIdx[i], 3 + j) = _PTPOptimizer->_noptJointCP(i, j);
			}
		}

		_qSpline = BSpline<-1, -1, -1>(_PTPOptimizer->_knot, cp);
		_qdotSpline = _qSpline.derivative();
		_qddotSpline = _qdotSpline.derivative();
	}

	void sharedResource::update(const VectorX & params)
	{
		//if (_params.size() == 0 || !_params.isApprox(params))
		//{
		//	_params = params;

		//	_tau.resize(_PTPOptimizer->_numOfGQSample, VectorX());
		//	_dtaudp.resize(_PTPOptimizer->_numOfGQSample, MatrixX());

		//	makeBSpline(_params);
		//	for (unsigned int i = 0; i < _PTPOptimizer->_numOfGQSample; i++)
		//	{
		//		_state[i]->setJointStatePos(_qSpline(_PTPOptimizer->GQ.getQueryPoints()[i]));
		//		_state[i]->setJointStateVel(_qdotSpline(_PTPOptimizer->GQ.getQueryPoints()[i]));
		//		_state[i]->setJointStateAcc(_qddotSpline(_PTPOptimizer->GQ.getQueryPoints()[i]));
		//		 
		//		_PTPOptimizer->_soc->solveInverseDynamics(*_state[i]);
		//		_tau[i] = _state[i]->getJointStateTorque();
		//		_dtaudp[i] = _PTPOptimizer->_soc->solveDiffInverseDynamics(*_state[i], _dqdp[i], _dqdotdp[i], _dqddotdp[i]);
		//	}
		//}
		if (_params.size() == 0 || !_params.isApprox(params))
		{
			_params = params;

			_tau.resize(_PTPOptimizer->_numOfData, VectorX());
			_dtaudp.resize(_PTPOptimizer->_numOfData, MatrixX());

			makeBSpline(_params);
			for (unsigned int i = 0; i < _PTPOptimizer->_numOfData; i++)
			{
				_state[i]->setJointStatePos(_qSpline(_PTPOptimizer->_tspan(i)));
				_state[i]->setJointStateVel(_qdotSpline(_PTPOptimizer->_tspan(i)));
				_state[i]->setJointStateAcc(_qddotSpline(_PTPOptimizer->_tspan(i)));

				_PTPOptimizer->_soc->solveInverseDynamics(*_state[i]);
				_tau[i] = _state[i]->getJointStateTorque();
				_dtaudp[i] = _PTPOptimizer->_soc->solveDiffInverseDynamics(*_state[i], _dqdp[i], _dqdotdp[i], _dqddotdp[i]);
			}
		}
	}

	void sharedResource::updatetau(const VectorX& params)
	{
		if (_params.size() == 0 || !_params.isApprox(params))
		{
			_params = params;

			_tau.resize(_PTPOptimizer->_numOfData, VectorX());

			makeBSpline(_params);
			for (unsigned int i = 0; i < _PTPOptimizer->_numOfData; i++)
			{
				_state[i]->setJointStatePos(_qSpline(_PTPOptimizer->_tspan(i)));
				_state[i]->setJointStateVel(_qdotSpline(_PTPOptimizer->_tspan(i)));
				_state[i]->setJointStateAcc(_qddotSpline(_PTPOptimizer->_tspan(i)));

				_PTPOptimizer->_soc->solveInverseDynamics(*_state[i]);
				_tau[i] = _state[i]->getJointStateTorque();
			}
		}
	}

	void sharedResource::updatedtaudp(const VectorX& params)
	{
		if (_params.size() == 0 || !_params.isApprox(params))
		{
			_params = params;

			_tau.resize(_PTPOptimizer->_numOfData, VectorX());
			_dtaudp.resize(_PTPOptimizer->_numOfData, MatrixX());

			makeBSpline(_params);
			for (unsigned int i = 0; i < _PTPOptimizer->_numOfData; i++)
			{
				_state[i]->setJointStatePos(_qSpline(_PTPOptimizer->_tspan(i)));
				_state[i]->setJointStateVel(_qdotSpline(_PTPOptimizer->_tspan(i)));
				_state[i]->setJointStateAcc(_qddotSpline(_PTPOptimizer->_tspan(i)));

				_PTPOptimizer->_soc->solveInverseDynamics(*_state[i]);
				_tau[i] = _state[i]->getJointStateTorque();
				_dtaudp[i] = _PTPOptimizer->_soc->solveDiffInverseDynamics(*_state[i], _dqdp[i], _dqdotdp[i], _dqddotdp[i]);
			}
		}
	}

	const vector<VectorX>& sharedResource::gettau(const VectorX & params)
	{
		//update(params);
		updatetau(params);
		return _tau;
	}

	const vector<MatrixX>& sharedResource::getdtaudp(const VectorX & params)
	{
		//update(params);
		updatedtaudp(params);
		return _dtaudp;
	}

	VectorX effortFunction::func(const VectorX & params) const
	{
		const vector<VectorX>& tau = _PTPOptimizer->_shared->gettau(params);
		const VectorX& weight = _PTPOptimizer->GQ.getWeights();
		VectorX fval = VectorX::Zero(1);

		for (unsigned int i = 0; i < _PTPOptimizer->_numOfGQSample; i++)
		{
			fval(0) += weight(i) * tau[i].squaredNorm();
		}
		return fval;
	}

	MatrixX effortFunction::Jacobian(const VectorX & params) const
	{
		const vector<VectorX>& tau = _PTPOptimizer->_shared->gettau(params);
		const std::vector<MatrixX>& dtaudp = _PTPOptimizer->_shared->getdtaudp(params);
		const VectorX& weight = _PTPOptimizer->GQ.getWeights();
		MatrixX jacobian = MatrixX::Zero(1, params.size());
		for (unsigned int i = 0; i < _PTPOptimizer->_numOfGQSample; i++)
		{
			jacobian += weight(i) * 2 * tau[i].transpose()*dtaudp[i];
		}
		return jacobian;
	}

	VectorX energyLossFunction::func(const VectorX& params) const
	{
		//const vector<VectorX>& tau = _PTPOptimizer->_shared->gettau(params);
		//const VectorX& weight = _PTPOptimizer->GQ.getWeights();

		//MotorJointPtr joint;
		//Real voltage, current;
		//VectorX fval = VectorX::Zero(1);
		//for (unsigned int j = 0; j < _PTPOptimizer->_soc->getNumOfJoint(); j++)
		//{
		//	for (unsigned int k = 0; k < _PTPOptimizer->_numOfGQSample; k++)
		//	{
		//		joint = _PTPOptimizer->_soc->getMotorJointPtr(j);

		//		current = 1.0 / (joint->getMotorConstant() * joint->getGearRatio()) * tau[k](j) +
		//			joint->getRotorInertia() * joint->getGearRatio() / joint->getMotorConstant() * _PTPOptimizer->_shared->getstate()[k]->getJointStateAcc(j);
		//		voltage = current * joint->getResistance() + joint->getBackEMFConstant() * joint->getGearRatio() * _PTPOptimizer->_shared->getstate()[k]->getJointStateVel(j);

		//		fval(0) += weight[k] * max(current * voltage, 0.0);
		//	}
		//}
		//return fval;

		const vector<VectorX>& tau = _PTPOptimizer->_shared->gettau(params);

		MotorJointPtr joint;
		Real voltage, current;
		VectorX fval = VectorX::Zero(1);
		for (unsigned int j = 0; j < _PTPOptimizer->_soc->getNumOfJoint(); j++)
		{
			for (unsigned int k = 0; k < _PTPOptimizer->_numOfData; k++)
			{
				joint = _PTPOptimizer->_soc->getMotorJointPtr(j);

				current = 1.0 / (joint->getMotorConstant() * joint->getGearRatio()) * tau[k](j) +
					joint->getRotorInertia() * joint->getGearRatio() * joint->getGearRatio() / joint->getMotorConstant() * _PTPOptimizer->_shared->getstate()[k]->getJointStateAcc(j);
				voltage = current * joint->getResistance() + joint->getBackEMFConstant() * joint->getGearRatio() * _PTPOptimizer->_shared->getstate()[k]->getJointStateVel(j);

				fval(0) += max(current * voltage, 0.0);
			}
		}
		fval(0) *= _PTPOptimizer->_dt;
		cout << "Energyloss val : " << fval(0) << endl;
		return fval;
	}

	//MatrixX energyLossFunction::Jacobian(const VectorX& params) const
	//{
	//	//const vector<VectorX>& tau = _PTPOptimizer->_shared->gettau(params);
	//	//const std::vector<MatrixX>& dtaudp = _PTPOptimizer->_shared->getdtaudp(params);
	//	//const VectorX& weight = _PTPOptimizer->GQ.getWeights();
	//	//MotorJointPtr joint;
	//	//Real voltage, current;
	//	//MatrixX jacobian = MatrixX::Zero(1, params.size());
	//	//for (unsigned int j = 0; j < _PTPOptimizer->_soc->getNumOfJoint(); j++)
	//	//{
	//	//	for (unsigned int k = 0; k < _PTPOptimizer->_numOfGQSample; k++)
	//	//	{
	//	//		joint = _PTPOptimizer->_soc->getMotorJointPtr(j);

	//	//		current = 1.0 / (joint->getMotorConstant() * joint->getGearRatio()) * tau[k](j) +
	//	//			joint->getRotorInertia() * joint->getGearRatio() / joint->getMotorConstant() * _PTPOptimizer->_shared->getstate()[k]->getJointStateAcc(j);
	//	//		voltage = current * joint->getResistance() + joint->getBackEMFConstant() * joint->getGearRatio() * _PTPOptimizer->_shared->getstate()[k]->getJointStateVel(j);

	//	//		if (RealBigger(current * voltage, 0.0))
	//	//		{
	//	//			MatrixX dcurrentdp = 1.0 / (joint->getMotorConstant() * joint->getGearRatio()) * dtaudp[k].row(j) +
	//	//				joint->getRotorInertia() * joint->getGearRatio() / joint->getMotorConstant() * _PTPOptimizer->_shared->getdqddotdp()[k].row(j);
	//	//			jacobian += weight[k] * ((current * joint->getResistance() + voltage) * dcurrentdp + current * joint->getBackEMFConstant() * joint->getGearRatio() * _PTPOptimizer->_shared->getdqdotdp()[k].row(j));	
	//	//		}
	//	//	}
	//	//}
	//	//return jacobian;

	//	const vector<VectorX>& tau = _PTPOptimizer->_shared->gettau(params);
	//	const std::vector<MatrixX>& dtaudp = _PTPOptimizer->_shared->getdtaudp(params);
	//	MotorJointPtr joint;
	//	Real voltage, current;
	//	MatrixX jacobian = MatrixX::Zero(1, params.size());
	//	for (unsigned int j = 0; j < _PTPOptimizer->_soc->getNumOfJoint(); j++)
	//	{
	//		for (unsigned int k = 0; k < _PTPOptimizer->_numOfData; k++)
	//		{
	//			joint = _PTPOptimizer->_soc->getMotorJointPtr(j);

	//			current = 1.0 / (joint->getMotorConstant() * joint->getGearRatio()) * tau[k](j) +
	//				joint->getRotorInertia() * joint->getGearRatio() / joint->getMotorConstant() * _PTPOptimizer->_shared->getstate()[k]->getJointStateAcc(j);
	//			voltage = current * joint->getResistance() + joint->getBackEMFConstant() * joint->getGearRatio() * _PTPOptimizer->_shared->getstate()[k]->getJointStateVel(j);

	//			if (RealBigger(current * voltage, 0.0))
	//			{
	//				MatrixX dcurrentdp = 1.0 / (joint->getMotorConstant() * joint->getGearRatio()) * dtaudp[k].row(j) +
	//					joint->getRotorInertia() * joint->getGearRatio() / joint->getMotorConstant() * _PTPOptimizer->_shared->getdqddotdp()[k].row(j);
	//				jacobian += ((current * joint->getResistance() + voltage) * dcurrentdp + current * joint->getBackEMFConstant() * joint->getGearRatio() * _PTPOptimizer->_shared->getdqdotdp()[k].row(j));
	//			}
	//		}
	//	}
	//	jacobian *= _PTPOptimizer->_dt;
	//	return jacobian;
	//}

	VectorX accelerationFunction::func(const VectorX& params) const
	{
		const vector<VectorX>& tau = _PTPOptimizer->_shared->gettau(params);
		const VectorX& querypoint = _PTPOptimizer->GQ.getQueryPoints();
		MatrixX qq(_PTPOptimizer->_soc->getNumOfJoint(), querypoint.size());

		const VectorX& weight = _PTPOptimizer->GQ.getWeights();
		VectorX fval = VectorX::Zero(1);

		for (unsigned int j = 0; j < _PTPOptimizer->_soc->getNumOfJoint(); j++)
		{
			for (unsigned int k = 0; k < _PTPOptimizer->_numOfGQSample; k++)
				fval(0) += weight[k] * _PTPOptimizer->_shared->_qddotSpline(querypoint[k]).squaredNorm();
		}
		return fval;
	}


	MatrixX accelerationFunction::Jacobian(const VectorX& params) const
	{
		const vector<VectorX>& tau = _PTPOptimizer->_shared->gettau(params);
		const VectorX& querypoint = _PTPOptimizer->GQ.getQueryPoints();
		const VectorX& weight = _PTPOptimizer->GQ.getWeights();
		MatrixX jacobian = MatrixX::Zero(1, params.size());

		VectorX v(_PTPOptimizer->_numOfOptCP);
		MatrixX dRdP(_PTPOptimizer->_shared->_dRdP.rows() + 2, _PTPOptimizer->_numOfOptCP);
		dRdP.setZero();
		dRdP.block(1, 0, _PTPOptimizer->_shared->_dRdP.rows(), _PTPOptimizer->_numOfOptCP) = _PTPOptimizer->_shared->_dRdP;

		MatrixX M(_PTPOptimizer->_soc->getNumOfJoint(), params.size()); M.setZero();

		std::vector<BSpline<-1, - 1, -1>> qs;
		for (int i = 0; i < v.size(); i++)
		{
			BSpline<-1, -1, -1> q(_PTPOptimizer->_shared->_qddotSpline.getKnots(), dRdP.col(i).transpose());
			qs.push_back(q);
		}

		for (unsigned int i = 0; i < _PTPOptimizer->_numOfGQSample; i++)
		{
			
			for (int j = 0; j < v.size(); j++)
				v(j) = qs[j](querypoint[i])(0);
			for (unsigned int j = 0; j < _PTPOptimizer->_numOfOptJoint; j++)
			{
				M.block(j, v.size()*j, 1, v.size()) = v.transpose();
			}
			jacobian += weight[i] * 2 * _PTPOptimizer->_shared->_qddotSpline(querypoint[i]).transpose() * M;
		}

		return jacobian;
	}


	VectorX NonlinearInequalityConstraint::func(const VectorX& params) const
	{
		const vector<VectorX>& tau = _PTPOptimizer->_shared->gettau(params);
		VectorX fval(_PTPOptimizer->_soc->getNumOfJoint() * 2);
		Real upper, lower;
		for (unsigned int i = 0; i < _PTPOptimizer->_soc->getNumOfJoint(); i++)
		{
			upper = RealMin;
			lower = RealMax;
			//for (unsigned int j = 0; j < _PTPOptimizer->_numOfGQSample; j++)
			for (unsigned int j = 0; j < _PTPOptimizer->_numOfData; j++)
			{
				if (upper < tau[j](i))
				{
					upper = tau[j](i);
				}
				if (lower > tau[j](i))
				{
					lower = tau[j](i);
				}
			}

			fval(i) = upper - _PTPOptimizer->_soc->getMotorJointPtr(i)->getLimitTorqueUpper();
			fval(_PTPOptimizer->_soc->getNumOfJoint() + i) = lower * (-1) + _PTPOptimizer->_soc->getMotorJointPtr(i)->getLimitTorqueLower();
		}
		return fval;
	}

	//MatrixX NonlinearInequalityConstraint::Jacobian(const VectorX& params) const
	//{
	//	const vector<VectorX>& tau = _PTPOptimizer->_shared->gettau(params);
	//	const std::vector<MatrixX>& dtaudp = _PTPOptimizer->_shared->getdtaudp(params);
	//	MatrixX jacobian(_PTPOptimizer->_soc->getNumOfJoint() * 2, params.size());
	//	Real upper, lower;
	//	unsigned int upperIdx, lowerIdx;
	//	for (unsigned int i = 0; i < _PTPOptimizer->_soc->getNumOfJoint(); i++)
	//	{
	//		upper = RealMin;
	//		lower = RealMax;
	//		//for (unsigned int j = 0; j < _PTPOptimizer->_numOfGQSample; j++)
	//		for (unsigned int j = 0; j < _PTPOptimizer->_numOfData; j++)
	//		{
	//			if (upper < tau[j](i))
	//			{
	//				upper = tau[j](i);
	//				upperIdx = j;
	//			}
	//			if (lower > tau[j](i))
	//			{
	//				lower = tau[j](i);
	//				lowerIdx = j;
	//			}
	//		}
	//		jacobian.row(i) = dtaudp[upperIdx].row(i);
	//		jacobian.row(_PTPOptimizer->_soc->getNumOfJoint() + i) = -dtaudp[lowerIdx].row(i);
	//	}
	//	return jacobian;
	//}
}
