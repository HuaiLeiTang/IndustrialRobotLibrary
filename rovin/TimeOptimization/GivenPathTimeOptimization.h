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
		std::vector<SE3, Eigen::aligned_allocator<SE3>> _givenPath;

		///< Minimum time
		Real _tm;

		std::vector<Real> _s; ///< trajectory parameter
		std::vector<Real> _sdot; ///< trajectory parameter's time derivative
		std::vector<Real> _sddot; ///< trajectory parameter's time 2nd derivative
		
		std::vector<VectorX, Eigen::aligned_allocator<VectorX>> _q; ///< joint angle
		std::vector<VectorX, Eigen::aligned_allocator<VectorX>> _qs; ///< joint angle's parameter derivative
		std::vector<VectorX, Eigen::aligned_allocator<VectorX>> _qss; ///< joint angle's paramter 2nd derivative

		SerialOpenChainPtr _socRobotPtr;

	public:
		GivenPathTimeOptimization();
		GivenPathTimeOptimization(const SerialOpenChainPtr& socRobotPtr, const std::vector<SE3, Eigen::aligned_allocator<SE3>>& givenPath);
		~GivenPathTimeOptimization();

		VectorX getq(const unsigned int i) const;
		VectorX getqs(const unsigned int i) const;
		VectorX getqss(const unsigned int i) const;
		SerialOpenChainPtr getSOCPtr();

		/*!
		* \brief Input GivenPath function
		*        Input path trajectory
		*		 change _givenPath, _q, _qs, _qss
		*/
		void InputGivenPath(const std::vector<SE3, Eigen::aligned_allocator<SE3>>& givenPath);

		void solveMinimumTimeOptimization();
		
		void solveInvKinAll();

		/*!
		* \brief Minium maximum acceleration finding function
		*        Input: parameter s(index), paramter's time derivative sdot
		*		 Output: minimum acceleration(min), maximum acceleration(max)
		*/
		void MinMaxAcc(const unsigned int index, const Real& sdot, /*output*/ Real& min, /*output*/ Real& max);


	};
}