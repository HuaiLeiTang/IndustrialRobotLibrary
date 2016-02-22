#include "PTPOptimizer.h"

namespace rovin {

	PTPoptimizer::PTPoptimizer(const SerialOpenChainPtr & soc, const std::vector<bool>& optJoint, 
		const unsigned int orderOfBSpline, const unsigned int numOfOptCP, const unsigned int numOfGQSample,
		const std::vector<Real>& tf, const std::vector<StatePtr>& constraintState)
	{
		LOGIF(constraintState.size() > 1, "PTPWayPointOptimizer::PTPWayPointOptimizer error : lack of number of constraint states.");
		LOGIF(constraintState.size() == (tf.size() + 1), "PTPWayPointOptimizer::PTPWayPointOptimizer error : tf size or constraint state size is wrong.");

		_soc = soc;
		_optJoint = optJoint;
		_orderOfBSpline = orderOfBSpline;
		_numOfOptCP = numOfOptCP;
		_numOfGQSample = numOfGQSample;
		_tf = tf;
		_constraintState = constraintState;

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




}
