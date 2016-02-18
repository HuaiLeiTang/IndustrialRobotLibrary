#include "PTPWayPointOptimization.h"

#include <memory>

using namespace std;

namespace rovin {

	PTPWayPointOptimizer::PTPWayPointOptimizer(const SerialOpenChainPtr & soc, const std::vector<bool>& optJoint, const unsigned int orderOfBSpline, const unsigned int numOfCP, const unsigned int numOfGQSample, const std::vector<Real>& tf, const std::vector<StatePtr>& constraintState)
	{
		LOGIF(constraintState.size() > 1, "PTPWayPointOptimizer::PTPWayPointOptimizer error : lack of number of constraint states.");
		LOGIF(constraintState.size() == (tf.size() + 1), "PTPWayPointOptimizer::PTPWayPointOptimizer error : tf size or constraint state size is wrong.");

		_numOfBspline = tf.size();
		_bsplineResult.resize(_numOfBspline);

		_wayPointOptimizer = PTPWayPointOptimizationPtr(new PTPWayPointOptimization(soc, optJoint, orderOfBSpline,
			numOfCP, numOfGQSample, tf[0], constraintState[0], constraintState[1], PTPWayPointOptimization::INITIAL_ALL_FINAL_POS));
	}

	void PTPWayPointOptimizer::generateTrajectory()
	{
		for (int i = 0; i < _numOfBspline; i++)
		{

			if (i == 0)
			{
				// first b-spline



			}
			else if (i == (_numOfBspline - 1))
			{
				// middle b-splines

			}
			else
			{
				// last b-spline


			}
		}
	}

	PTPWayPointOptimization::PTPWayPointOptimization(const SerialOpenChainPtr & soc, const vector<bool>& optJoint,
		const unsigned int orderOfBSpline, const unsigned int numOfOptCP, const unsigned int numOfGQSample,
		const Real tf, const StatePtr& initialState, const StatePtr& finalState, ConstraintCondition constraintCondition)
	{
		_soc = soc;
		_optJoint = optJoint;
		_orderOfBSpline = orderOfBSpline;
		_numOfOptCP = numOfOptCP;
		_numOfGQSample = numOfGQSample;
		_tf = tf;
		_initialState = initialState;
		_finalState = finalState;
		_constraintCondition = constraintCondition;


		switch (constraintCondition)
		{
		case INITIAL_ALL_FINAL_ALL:
			_numOfConstraintCP = 6;
			break;

		case INITIAL_ALL_FINAL_POS:
			_numOfConstraintCP = 4;
			break;

		default:
			break;
		}

		// insert optimizing joint index and non-optimizing joint index
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

	void PTPWayPointOptimization::setInitialState(const StatePtr & initState)
	{
		_initialState = initState;
	}

	void PTPWayPointOptimization::setFinalState(const StatePtr & finalState)
	{
		_finalState = finalState;
	}

	void PTPWayPointOptimization::setConstraintCondition(const ConstraintCondition constraintCondition)
	{
		_constraintCondition = constraintCondition;
	}

	void PTPWayPointOptimization::setFinalTime(const Real tf)
	{
		_tf = tf;
	}

	void PTPWayPointOptimization::makeBSplineKnot()
	{
		LOGIF(_orderOfBSpline + _numOfOptCP + _numOfConstraintCP > 2 * _orderOfBSpline, "The number of control points is not enough.");

		//_knot.resize(_orderOfBSpline + _numOfOptCP + _numOfConstraintCP);
		//for (unsigned int i = 0; i < _orderOfBSpline; i++)
		//{
		//	_knot[i] = 0;
		//	_knot[_knot.size() - i - 1] = _tf;
		//}

		//Real delta = _tf / (_numOfOptCP + _numOfConstraintCP - _orderOfBSpline + 1);
		//for (unsigned int i = 0; i < _numOfOptCP + _numOfConstraintCP - _orderOfBSpline; i++)
		//{
		//	_knot[_orderOfBSpline + i] = delta * (i + 1);
		//}

		_knot.resize(_orderOfBSpline + _numOfOptCP + _numOfConstraintCP + 1);
		for (unsigned int i = 0; i < _orderOfBSpline + 1; i++)
		{
			_knot[i] = 0;
			_knot[_knot.size() - i - 1] = _tf;
		}

		Real delta = _tf / (_numOfOptCP + _numOfConstraintCP - _orderOfBSpline);
		for (unsigned int i = 0; i < _numOfOptCP + _numOfConstraintCP - _orderOfBSpline; i++)
		{
			_knot[_orderOfBSpline + 1 + i] = delta * (i + 1);
		}
	}

	void PTPWayPointOptimization::makeBoundaryCondition()
	{
		Real delta = _tf / (_numOfOptCP + _numOfConstraintCP - _orderOfBSpline);

		if(_constraintCondition == INITIAL_ALL_FINAL_ALL)
		{ 
			_initialCP.resize(3);
			_initialCP[0] = _initialState->getJointStatePos();
			_initialCP[1] = delta / (_orderOfBSpline - 1)*_initialState->getJointStateVel() + _initialCP[0];
			_initialCP[2] = 2 * delta*(delta / (_orderOfBSpline - 1) / (_orderOfBSpline - 2)*_initialState->getJointStateAcc() + _initialCP[1] * (1 / (2 * delta) + 1 / delta) - _initialCP[0] / delta);

			_finalCP.resize(3);
			_finalCP[0] = _finalState->getJointStatePos();
			_finalCP[1] = -delta / (_orderOfBSpline - 1)*_finalState->getJointStateVel() + _finalCP[0];
			_finalCP[2] = 2 * delta*(delta / (_orderOfBSpline - 1) / (_orderOfBSpline - 2)*_finalState->getJointStateAcc() + _finalCP[1] * (1 / (2 * delta) + 1 / delta) - _finalCP[0] / delta);
		}
		else if(_constraintCondition == INITIAL_ALL_FINAL_POS)
		{
			_initialCP.resize(3);
			_initialCP[0] = _initialState->getJointStatePos();
			_initialCP[1] = delta / (_orderOfBSpline - 1)*_initialState->getJointStateVel() + _initialCP[0];
			_initialCP[2] = 2 * delta*(delta / (_orderOfBSpline - 1) / (_orderOfBSpline - 2)*_initialState->getJointStateAcc() + _initialCP[1] * (1 / (2 * delta) + 1 / delta) - _initialCP[0] / delta);

			_finalCP.resize(1);
			_finalCP[0] = _finalState->getJointStatePos();
		}
		
	}


	const BSpline<-1, -1, -1>& PTPWayPointOptimization::getResultBSpline() const
	{
		return _shared->_qSpline;
	}

	void PTPWayPointOptimization::makeNonOptJointCP()
	{
		_noptJointCP.resize(_soc->getNumOfJoint() - _numOfOptJoint, _numOfOptCP);
		for (unsigned int i = 0; i < _soc->getNumOfJoint() - _numOfOptJoint; i++)
		{
			for (unsigned int j = 0; j < _numOfOptCP; j++)
			{
				if (_constraintCondition == INITIAL_ALL_FINAL_ALL)
				{
					_noptJointCP(i, j) = (_finalCP[2](_noptJointIdx[i]) - _initialCP[2](_noptJointIdx[i])) / (_numOfOptCP + 1) * (j + 1) + _initialCP[2](_noptJointIdx[i]);
				}
				else if (_constraintCondition == INITIAL_ALL_FINAL_POS)
				{
					_noptJointCP(i, j) = (_finalCP[0](_noptJointIdx[i]) - _initialCP[2](_noptJointIdx[i])) / (_numOfOptCP + 1) * (j + 1) + _initialCP[2](_noptJointIdx[i]);
				}
			}
		}
	}

	void PTPWayPointOptimization::makeObjectiveFunction()
	{
		_objectFunc = FunctionPtr(new effortFunctionWayPoint(this));
	}

	/////////////////////


	void PTPWayPointOptimization::makeIneqConstraintFunction()
	{
		_nonlinearIneqFunc = FunctionPtr(new NonlinearInequalityConstraintWayPoint(this));

		// make linear inequality function
		// 곱하기 2 : upper bound & lower bound
		unsigned int np = _shared->_dPdP.rows();
		unsigned int nq = _shared->_dQdP.rows();
		unsigned int nr = _shared->_dRdP.rows();
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

			b.block((np + nq + nr) * 2 * i + np * 2, 0, nq, 1) = _shared->_Q[i] - VectorX::Ones(nq)*_soc->getMotorJointPtr(_optJointIdx[i])->getLimitVelUpper();
			b.block((np + nq + nr) * 2 * i + np * 2 + nq, 0, nq, 1) = -_shared->_Q[i] + VectorX::Ones(nq)*_soc->getMotorJointPtr(_optJointIdx[i])->getLimitVelLower();

			b.block((np + nq + nr) * 2 * i + (np + nq) * 2, 0, nr, 1) = _shared->_R[i] - VectorX::Ones(nr)*_soc->getMotorJointPtr(_optJointIdx[i])->getLimitAccUpper();
			b.block((np + nq + nr) * 2 * i + (np + nq) * 2 + nr, 0, nr, 1) = -_shared->_R[i] + VectorX::Ones(nr)*_soc->getMotorJointPtr(_optJointIdx[i])->getLimitAccLower();
		}
		_linearIneqFunc = FunctionPtr(new AffineFunction(A, b));


		_IneqFunc = FunctionPtr(new AugmentedFunction());
		static_pointer_cast<AugmentedFunction>(_IneqFunc)->addFunction(_nonlinearIneqFunc);
		static_pointer_cast<AugmentedFunction>(_IneqFunc)->addFunction(_linearIneqFunc);
	}

	void PTPWayPointOptimization::generateTrajectory()
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

		_shared = sharedResourceWayPointPtr(new sharedResourceWayPoint(this));
		LOG("Complete shared initialization.");
		//cout << _shared->_dPdP << endl;
		//cout << _shared->_dQdP << endl;
		//cout << _shared->_dRdP << endl;

		makeNonOptJointCP();

		makeObjectiveFunction();
		makeIneqConstraintFunction();
		LOG("Optimization ready.");
		//cout << "Nopt Joint CP: " << endl;
		//cout << _noptJointCP << endl;
		//cout << "LinearInequalityCondition(A): " << endl << static_pointer_cast<AffineFunction>(_linearIneqFunc)->getA() << endl;
		//cout << "LinearInequalityCondition(b): " << endl << static_pointer_cast<AffineFunction>(_linearIneqFunc)->getb() << endl;

		// Initial Guess (initial control points)
		VectorX initX(_numOfOptJoint * _numOfOptCP);
		for (unsigned int i = 0; i < _numOfOptJoint; i++)
		{
			for (unsigned int j = 0; j < _numOfOptCP; j++)
			{
				if (_constraintCondition == INITIAL_ALL_FINAL_ALL)
				{
					initX(_numOfOptCP * i + j) = (_finalCP[2](_optJointIdx[i]) - _initialCP[2](_optJointIdx[i])) / (_numOfOptCP + 1) * (j + 1) + _initialCP[2](_optJointIdx[i]);
				}
				else if (_constraintCondition == INITIAL_ALL_FINAL_POS)
				{
					initX(_numOfOptCP * i + j) = (_finalCP[0](_optJointIdx[i]) - _initialCP[2](_optJointIdx[i])) / (_numOfOptCP + 1) * (j + 1) + _initialCP[2](_optJointIdx[i]);
				}
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

		_optimizer.setObjectiveFunction(_objectFunc);
		_optimizer.setInequalityConstraint(_IneqFunc);
		LOG("Start optimization.");
		_optimizer.solve(initX);
		LOG("Finish optimization.");

		_shared->makeBSpline(_optimizer.resultX);

		//cout << "X" << endl;
		//cout << _optimizer.resultX << endl;
		//cout << "knot" << endl;
		//cout << _knot << endl;
		//cout << "controlpoint" << endl;
		//cout << _shared->_qSpline.getControlPoints() << endl;
		//cout << "Inequality : " << _IneqFunc->func(_optimizer.resultX) << endl;
		//cout << "f(X) : " << _optimizer.resultFunc << endl;
	}

	// class sharedResource

	sharedResourceWayPoint::sharedResourceWayPoint(PTPWayPointOptimization* PTPOptimizer)
	{
		_PTPOptimizer = PTPOptimizer;
		_state.resize(_PTPOptimizer->_numOfGQSample);
		for (unsigned int i = 0; i < _state.size(); i++)
		{
			_state[i] = _PTPOptimizer->_soc->makeState();
		}

		// Calculate _dqdp,dqdotdp, dqddotdq & Calculate _dPdP,_dQdP, _dRdP
		MatrixX cp(1, _PTPOptimizer->_numOfOptCP + _PTPOptimizer->_numOfConstraintCP); /// control points including boundary control points
		bool checkMatrixSize = false;
		_dqdp.resize(_PTPOptimizer->_numOfGQSample, MatrixX::Zero(_PTPOptimizer->_soc->getNumOfJoint(), _PTPOptimizer->_numOfOptCP * _PTPOptimizer->_numOfOptJoint));
		_dqdotdp.resize(_PTPOptimizer->_numOfGQSample, MatrixX::Zero(_PTPOptimizer->_soc->getNumOfJoint(), _PTPOptimizer->_numOfOptCP * _PTPOptimizer->_numOfOptJoint));
		_dqddotdp.resize(_PTPOptimizer->_numOfGQSample, MatrixX::Zero(_PTPOptimizer->_soc->getNumOfJoint(), _PTPOptimizer->_numOfOptCP * _PTPOptimizer->_numOfOptJoint));
		for (unsigned int i = 0; i < _PTPOptimizer->_numOfOptCP; i++)
		{
			cp.setZero();
			cp(0, 3 + i) = 1.0; // boundary control point 3개 이후에 1 값을 넣어준다
			_qSpline = BSpline<-1, -1, -1>(_PTPOptimizer->_knot, cp);
			_qdotSpline = _qSpline.derivative();
			_qddotSpline = _qdotSpline.derivative();

			if (!checkMatrixSize)
			{
				if (_PTPOptimizer->_constraintCondition == PTPWayPointOptimization::INITIAL_ALL_FINAL_ALL)
				{
					_dPdP.resize(_qSpline.getControlPoints().cols() - 6, _PTPOptimizer->_numOfOptCP);
					_dQdP.resize(_qdotSpline.getControlPoints().cols() - 4, _PTPOptimizer->_numOfOptCP);
					_dRdP.resize(_qddotSpline.getControlPoints().cols() - 2, _PTPOptimizer->_numOfOptCP);
				}
				else if (_PTPOptimizer->_constraintCondition == PTPWayPointOptimization::INITIAL_ALL_FINAL_POS)
				{
					_dPdP.resize(_qSpline.getControlPoints().cols() - 4, _PTPOptimizer->_numOfOptCP);
					_dQdP.resize(_qdotSpline.getControlPoints().cols() - 2, _PTPOptimizer->_numOfOptCP);
					_dRdP.resize(_qddotSpline.getControlPoints().cols() - 1, _PTPOptimizer->_numOfOptCP);
				}
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
			// control point에 initial boundary 값 3개 넣어주고
			cp(0) = _PTPOptimizer->_initialCP[0](_PTPOptimizer->_optJointIdx[i]);
			cp(1) = _PTPOptimizer->_initialCP[1](_PTPOptimizer->_optJointIdx[i]);
			cp(2) = _PTPOptimizer->_initialCP[2](_PTPOptimizer->_optJointIdx[i]);
			// control point에 final boundary 값 3개 넣어주고
			if (_PTPOptimizer->_constraintCondition == PTPWayPointOptimization::INITIAL_ALL_FINAL_ALL)
			{
				cp(_PTPOptimizer->_numOfOptCP + 3) = _PTPOptimizer->_finalCP[2](_PTPOptimizer->_optJointIdx[i]);
				cp(_PTPOptimizer->_numOfOptCP + 4) = _PTPOptimizer->_finalCP[1](_PTPOptimizer->_optJointIdx[i]);
				cp(_PTPOptimizer->_numOfOptCP + 5) = _PTPOptimizer->_finalCP[0](_PTPOptimizer->_optJointIdx[i]);
			}
			else if (_PTPOptimizer->_constraintCondition == PTPWayPointOptimization::INITIAL_ALL_FINAL_POS)
			{
				cp(_PTPOptimizer->_numOfOptCP + 3) = _PTPOptimizer->_finalCP[0](_PTPOptimizer->_optJointIdx[i]);
			}

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

	void sharedResourceWayPoint::makeBSpline(const VectorX& params)
	{
		MatrixX cp(_PTPOptimizer->_soc->getNumOfJoint(), _PTPOptimizer->_numOfOptCP + _PTPOptimizer->_numOfConstraintCP);

		// Boundary
		if (_PTPOptimizer->_constraintCondition == PTPWayPointOptimization::INITIAL_ALL_FINAL_ALL)
		{
			cp.col(0) = _PTPOptimizer->_initialCP[0];
			cp.col(1) = _PTPOptimizer->_initialCP[1];
			cp.col(2) = _PTPOptimizer->_initialCP[2];
			cp.col(cp.cols() - 1) = _PTPOptimizer->_finalCP[0];
			cp.col(cp.cols() - 2) = _PTPOptimizer->_finalCP[1];
			cp.col(cp.cols() - 3) = _PTPOptimizer->_finalCP[2];
		}
		else if (_PTPOptimizer->_constraintCondition == PTPWayPointOptimization::INITIAL_ALL_FINAL_POS)
		{
			cp.col(0) = _PTPOptimizer->_initialCP[0];
			cp.col(1) = _PTPOptimizer->_initialCP[1];
			cp.col(2) = _PTPOptimizer->_initialCP[2];
			cp.col(cp.cols() - 1) = _PTPOptimizer->_finalCP[0];
		}
		
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

	void sharedResourceWayPoint::update(const VectorX & params)
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

	const vector<VectorX>& sharedResourceWayPoint::gettau(const VectorX & params)
	{
		update(params);
		return _tau;
	}

	const vector<MatrixX>& sharedResourceWayPoint::getdtaudp(const VectorX & params)
	{
		update(params);
		return _dtaudp;
	}

	///

	VectorX effortFunctionWayPoint::func(const VectorX & params) const
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

	MatrixX effortFunctionWayPoint::Jacobian(const VectorX & params) const
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

	VectorX NonlinearInequalityConstraintWayPoint::func(const VectorX& params) const
	{
		//
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

	MatrixX NonlinearInequalityConstraintWayPoint::Jacobian(const VectorX& params) const
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