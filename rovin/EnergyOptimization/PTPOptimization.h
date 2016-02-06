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

namespace rovin {

	class PTPOptimization;
	class effortFunction;
	class energyLossFunction;
	class EqualityConstraint;
	class InequalityConstraint;


	class PTPOptimization
	{
	public:
		enum ObjectiveFunctionType { effort, energyloss };
		

	private:
		NonlinearOptimization _optimizer;

	public:
		




	};

	// 
	class effortFunction : Function
	{
		effortFunction() {}

		VectorX func(const VectorX& x) const;
		MatrixX Jacobian(const VectorX& x) const;
	};




	class energyLossFunction : Function
	{
		energyLossFunction() {}

	};

	class equalityConstraint :Function
	{
		equalityConstraint() {}

	};

	class InequalityConstraint :Function
	{
		InequalityConstraint() {}

	};


}
