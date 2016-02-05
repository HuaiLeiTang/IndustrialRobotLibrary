/*!
*	\file	GivenPathTimeOptimization.h
*	\date	2016.02.05
*	\author	Wooyoung (wykim1989@gmail.com)
*	\brief	Given Path Time Optimization
*          this class has GivenPathTimeOptimization
*/

#pragma once

#include <memory>
#include <vector>

#include <rovin/Math/Common.h>
#include <rovin/Math/LieGroup.h>
#include <rovin/Dynamics/SerialOpenChain.h>

namespace rovin
{
	class GivenPathTimeOptimization
	{
	private:
		///< Given endeffector trajctory 
		std::vector<SE3> _givenPath;

		///< Minimum time
		Real _tm;

		std::vector<Real> _s; ///< trajectory parameter
		std::vector<Real> _sdot; ///< trajectory parameter's time derivative
		std::vector<Real> _sddot; ///< trajectory parameter's time 2nd derivative
		
		std::vector<VectorX> _q; ///< joint angle
		std::vector<VectorX> _qs; ///< joint angle's parameter derivative
		std::vector<VectorX> _qss; ///< joint angle's paramter 2nd derivative

		SerialOpenChainPtr _socRobotPtr;

	public:
		GivenPathTimeOptimization();
		GivenPathTimeOptimization(const std::vector<SE3>& givenPath);
		~GivenPathTimeOptimization();

		VectorX getq(const unsigned int i);
		VectorX getqs(const unsigned int i);
		VectorX getqss(const unsigned int i);

		void setGivenPath(const std::vector<SE3>& givenPath);
		void solveMinimumTimeOptimization();

		VectorX MinAcc(const Real& s, const Real& sdot);
		VectorX MaxAcc(const Real& s, const Real& sdot);


	};
}