#include "PTPOptimization.h"

#include <time.h>
#include <memory>

using namespace std;

namespace rovin{

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
		makeBSplineKnot();
		LOG("Complete making B-Spline knots.");
		//cout << "knot : " << _knot.transpose() << endl;
		makeBoundaryCondition();
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

		makeNonOptJointCP();
		makeObjectiveFunction();
		if (_optType == OptimizationType::nlopt)
		{
			makeIneqConstraintFunction_nlopt();
		}
		else if ( (_optType == OptimizationType::GCMMA) || (_optType == OptimizationType::GCMMA_TR) || (_optType == OptimizationType::GCMMA_GD))
		{
			//makeIneqConstraintFunction_nlopt();
			makeIneqConstraintFunction_MMA();
		}
		
		LOG("Optimization ready.");
		//cout << "Nopt Joint CP: " << endl;
		//cout << _noptJointCP << endl;
		//cout << "LinearInequalityCondition(A): " << endl << static_pointer_cast<AffineFunction>(_linearIneqFunc)->getA() << endl;
		//cout << "LinearInequalityCondition(b): " << endl << static_pointer_cast<AffineFunction>(_linearIneqFunc)->getb() << endl;

		// Initial Guess
		VectorX initX(_numOfOptJoint * _numOfOptCP);
		for (unsigned int i = 0; i < _numOfOptJoint; i++)
		{
			for (unsigned int j = 0; j < _numOfOptCP; j++)
			{
				initX(_numOfOptCP * i + j) = (_finalCP[2](_optJointIdx[i]) - _initialCP[2](_optJointIdx[i])) / (_numOfOptCP + 1) * (j + 1) + _initialCP[2](_optJointIdx[i]);
			}
		}

		//initX << 0.556274, 0.11917, -0.328995, -0.814881, 0.603234, 0.0932657, -0.130955, -0.064882, 1.14852, 0.774952, 0.2058, -0.233091;
		//initX << 0.419949, 0.275863, -0.263421, -0.847695, -0.218918, -0.226281, -0.224308, -0.177036, 1.22172, 1.20281, 0.169464, -0.362267;
		LOG("Initial guess ready.");
		//cout << "Initial X: " << endl;
		//cout << initX << endl;
		//cout << "f(X): " << endl;
		//cout << _objectFunc->func(initX) << endl;
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
		}
		else if ((_optType == OptimizationType::GCMMA) || (_optType == OptimizationType::GCMMA_TR) || (_optType == OptimizationType::GCMMA_GD))
		{
			VectorX minX(initX.size()), maxX(initX.size());
			//for (int iii = 0; iii < _numOfOptCP; iii++)
			//{
			//	for (int jjj = 0; jjj < _numOfOptJoint; jjj++)
			//	{
			//		minX(iii * _numOfOptJoint + jjj) = _soc->getMotorJointPtr(_optJointIdx[jjj])->getLimitPosLower();
			//		maxX(iii * _numOfOptJoint + jjj) = _soc->getMotorJointPtr(_optJointIdx[jjj])->getLimitPosUpper();
			//	}
			//}
			for (unsigned int iii = 0; iii < _numOfOptJoint; iii++)
			{
				for (unsigned int jjj = 0; jjj < _numOfOptCP; jjj++)
				{
					minX(iii * _numOfOptCP + jjj) = _soc->getMotorJointPtr(_optJointIdx[iii])->getLimitPosLower();
					maxX(iii * _numOfOptCP + jjj) = _soc->getMotorJointPtr(_optJointIdx[iii])->getLimitPosUpper();
				}
			}

			//cout << "minX" << endl << minX << endl << endl;
			//cout << "maxX" << endl << maxX << endl << endl;


			_GCMMAoptimizer.initialize(initX.size(), _IneqFunc->func(initX).size());
			_GCMMAoptimizer.setMinMax(minX, maxX);
			_GCMMAoptimizer.setObjectiveFunction(_objectFunc);
			_GCMMAoptimizer.setInequalityConstraint(_IneqFunc);
			LOG("Start optimization.");
			clock_t time = clock();
			//for (int i = 0; i < 20; i ++)
			//	_GCMMAoptimizer.solve(initX);
			if (_optType == OptimizationType::GCMMA)
			{
				//for (int l = 0; l < 100; l++)
					_GCMMAoptimizer.solve(initX);
			}
			else if (_optType == OptimizationType::GCMMA_TR)
			{
				for (int l = 0; l < 100; l++)
					_GCMMAoptimizer.TR_solve(initX);
			}
			else if (_optType == OptimizationType::GCMMA_GD)
			{
				//for (int l = 0; l < 100; l ++)
					_GCMMAoptimizer.GD_solve(initX);
			}
			LOG("Finish optimization.");
			cout << "------------------------------------" << endl;
			cout << "computation time : " << (clock() - time) << endl << endl;
			cout << "X : " << endl << _GCMMAoptimizer.resultX << endl << endl;
			cout << "control points" << endl << _shared->_qSpline.getControlPoints() << endl << endl;
			cout << "Value of objective function : " << _GCMMAoptimizer.resultFunc << endl << endl;
			//if (_optType == OptimizationType::GCMMA)
			//{
			//	cout << "_sublam : " << endl << _GCMMAoptimizer._sublam << endl << endl;
			//	cout << "_suby : " << endl << _GCMMAoptimizer._suby << endl << endl;
			//}
			//else if (_optType == OptimizationType::GCMMA_TR)
			//{
			//	cout << "_sublam : " << endl << _GCMMAoptimizer._TR_sublam << endl << endl;
			//	cout << "_suby : " << endl << _GCMMAoptimizer._TR_suby << endl << endl;
			//}
			//else if (_optType == OptimizationType::GCMMA_GD)
			//{
			//	cout << "_sublam : " << endl << _GCMMAoptimizer._GD_sublam << endl << endl;
			//	cout << "_suby : " << endl << _GCMMAoptimizer._TR_suby << endl << endl;
			//}
				
		}
	}
	
	// class sharedResource

	sharedResource::sharedResource(PTPOptimization* PTPOptimizer)
	{
		_PTPOptimizer = PTPOptimizer;
		_state.resize(_PTPOptimizer->_numOfGQSample);
		for (unsigned int i = 0; i < _state.size(); i++)
		{
			_state[i] = _PTPOptimizer->_soc->makeState();
		}

		// Calculate _dqdp,dqdotdp, dqddotdq
		// Calculate _dPdP,_dQdP, _dRdP
		MatrixX cp(1, _PTPOptimizer->_numOfOptCP + 6);
		bool checkMatrixSize = false;
		_dqdp.resize(_PTPOptimizer->_numOfGQSample, MatrixX::Zero(_PTPOptimizer->_soc->getNumOfJoint(), _PTPOptimizer->_numOfOptCP * _PTPOptimizer->_numOfOptJoint));
		_dqdotdp.resize(_PTPOptimizer->_numOfGQSample, MatrixX::Zero(_PTPOptimizer->_soc->getNumOfJoint(), _PTPOptimizer->_numOfOptCP * _PTPOptimizer->_numOfOptJoint));
		_dqddotdp.resize(_PTPOptimizer->_numOfGQSample, MatrixX::Zero(_PTPOptimizer->_soc->getNumOfJoint(), _PTPOptimizer->_numOfOptCP * _PTPOptimizer->_numOfOptJoint));
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

			for (unsigned int j = 0; j < _PTPOptimizer->_numOfGQSample; j++)
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
		//cout << "_dPdP" << endl << _dPdP << endl;
		//cout << "_dQdP" << endl << _dQdP << endl;
		//cout << "_dRdP" << endl << _dRdP << endl;

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
		if (_params.size() == 0 || !_params.isApprox(params))
		{
			_params = params;

			_tau.resize(_PTPOptimizer->_numOfGQSample, VectorX());
			_dtaudp.resize(_PTPOptimizer->_numOfGQSample, MatrixX());

			makeBSpline(_params);
			for (unsigned int i = 0; i < _PTPOptimizer->_numOfGQSample; i++)
			{
				_state[i]->setJointStatePos(_qSpline(_PTPOptimizer->GQ.getQueryPoints()[i]));
				_state[i]->setJointStateVel(_qdotSpline(_PTPOptimizer->GQ.getQueryPoints()[i]));
				_state[i]->setJointStateAcc(_qddotSpline(_PTPOptimizer->GQ.getQueryPoints()[i]));
				 
				_PTPOptimizer->_soc->solveInverseDynamics(*_state[i]);
				_tau[i] = _state[i]->getJointStateTorque();
				_dtaudp[i] = _PTPOptimizer->_soc->solveDiffInverseDynamics(*_state[i], _dqdp[i], _dqdotdp[i], _dqddotdp[i]);
			}
		}
	}

	const vector<VectorX>& sharedResource::gettau(const VectorX & params)
	{
		update(params);
		return _tau;
	}

	const vector<MatrixX>& sharedResource::getdtaudp(const VectorX & params)
	{
		update(params);
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
		const vector<VectorX>& tau = _PTPOptimizer->_shared->gettau(params);
		const VectorX& weight = _PTPOptimizer->GQ.getWeights();
		MotorJointPtr joint;
		Real voltage, current;
		VectorX fval = VectorX::Zero(1);
		for (unsigned int j = 0; j < _PTPOptimizer->_soc->getNumOfJoint(); j++)
		{
			for (unsigned int k = 0; k < _PTPOptimizer->_numOfGQSample; k++)
			{
				joint = _PTPOptimizer->_soc->getMotorJointPtr(j);

				current = 1.0 / (joint->getMotorConstant() * joint->getGearRatio()) * tau[k](j) +
					joint->getRotorInertia() * joint->getGearRatio() / joint->getMotorConstant() * _PTPOptimizer->_shared->getstate()[k]->getJointStateAcc(j);
				voltage = current * joint->getResistance() + joint->getBackEMFConstant() * joint->getGearRatio() * _PTPOptimizer->_shared->getstate()[k]->getJointStateVel(j);

				fval(0) += weight[k] * max(current * voltage, 0.0);
			}
		}
		return fval;
	}

	MatrixX energyLossFunction::Jacobian(const VectorX& params) const
	{
		const vector<VectorX>& tau = _PTPOptimizer->_shared->gettau(params);
		const std::vector<MatrixX>& dtaudp = _PTPOptimizer->_shared->getdtaudp(params);
		const VectorX& weight = _PTPOptimizer->GQ.getWeights();
		MotorJointPtr joint;
		Real voltage, current;
		MatrixX jacobian = MatrixX::Zero(1, params.size());
		for (unsigned int j = 0; j < _PTPOptimizer->_soc->getNumOfJoint(); j++)
		{
			for (unsigned int k = 0; k < _PTPOptimizer->_numOfGQSample; k++)
			{
				joint = _PTPOptimizer->_soc->getMotorJointPtr(j);

				current = 1.0 / (joint->getMotorConstant() * joint->getGearRatio()) * tau[k](j) +
					joint->getRotorInertia() * joint->getGearRatio() / joint->getMotorConstant() * _PTPOptimizer->_shared->getstate()[k]->getJointStateAcc(j);
				voltage = current * joint->getResistance() + joint->getBackEMFConstant() * joint->getGearRatio() * _PTPOptimizer->_shared->getstate()[k]->getJointStateVel(j);

				if (RealBigger(current * voltage, 0.0))
				{
					MatrixX dcurrentdp = 1.0 / (joint->getMotorConstant() * joint->getGearRatio()) * dtaudp[k].row(j) +
						joint->getRotorInertia() * joint->getGearRatio() / joint->getMotorConstant() * _PTPOptimizer->_shared->getdqddotdp()[k].row(j);
					jacobian += weight[k] * ((current * joint->getResistance() + voltage) * dcurrentdp + current * joint->getBackEMFConstant() * joint->getGearRatio() * _PTPOptimizer->_shared->getdqdotdp()[k].row(j));	
				}
			}
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
			for (unsigned int j = 0; j < _PTPOptimizer->_numOfGQSample; j++)
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

	MatrixX NonlinearInequalityConstraint::Jacobian(const VectorX& params) const
	{
		const vector<VectorX>& tau = _PTPOptimizer->_shared->gettau(params);
		const std::vector<MatrixX>& dtaudp = _PTPOptimizer->_shared->getdtaudp(params);
		MatrixX jacobian(_PTPOptimizer->_soc->getNumOfJoint() * 2, params.size());
		Real upper, lower;
		unsigned int upperIdx, lowerIdx;
		for (unsigned int i = 0; i < _PTPOptimizer->_soc->getNumOfJoint(); i++)
		{
			upper = RealMin;
			lower = RealMax;
			for (unsigned int j = 0; j < _PTPOptimizer->_numOfGQSample; j++)
			{
				if (upper < tau[j](i))
				{
					upper = tau[j](i);
					upperIdx = j;
				}
				if (lower > tau[j](i))
				{
					lower = tau[j](i);
					lowerIdx = j;
				}
			}

			jacobian.row(i) = dtaudp[upperIdx].row(i);
			jacobian.row(_PTPOptimizer->_soc->getNumOfJoint() + i) = -dtaudp[lowerIdx].row(i);
		}
		return jacobian;
	}
}
