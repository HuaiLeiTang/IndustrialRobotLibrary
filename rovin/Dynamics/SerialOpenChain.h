/*!
 *	\file	SerialOpenChain.h
 *	\date	2016.01.22
 *	\author	Youngsuk (crazyhys@gmail.com)
 *	\brief	SerialOpenChain class
 *          this class has robot properties including assembly
 *          Robot assembly is implemented in this class
*/

#pragma once

#include <vector>
#include <memory>

#include <rovin\Utils\Diagnostic.h>

#include "Link.h"
#include "MotorJoint.h"
#include "State.h"

namespace rovin
{
	class SerialOpenChain;
	class SerialOpenChainLink;
	class SerialOpenChainJoint;
	class Mate;
	typedef std::shared_ptr<SerialOpenChain> SerialOpenChainPtr;

	class SerialOpenChain
	{
	public:
		/// forward kinematics option
		enum JOINT_KINEMATICS_OPTION
		{
			//JOINT_TRANSFORM = 1 << 0, ///< option for forward kinematics, endeffector SE3 
			//JOINT_JACOBIAN = 1 << 1, ///< option for differenctial forward kinematics, Link velocity
			//JOINT_JACOBIANDOT = 1 << 2 ///< option for differenctial forward kinematics, Link acceleration
		};

	private:
		/*!
		* \brief SerialOpenChain class member variable
		*/
		std::vector<LinkPtr> _linkPtr;
		std::vector<MotorJointPtr> _motorJointPtr;
		std::vector<Mate> _Mate;

		std::vector<SerialOpenChainLink, Eigen::aligned_allocator< SerialOpenChainLink >> _socLink;
		std::vector<SerialOpenChainJoint, Eigen::aligned_allocator< SerialOpenChainLink >> _socJoint;

		bool _complete;

	public:
		/*!
		* \brief SerialOpenChain class member functions
		*/

		// constructor & destructor
		SerialOpenChain();
		~SerialOpenChain();
		
		// get-function
		LinkPtr getLinkPtr(const unsigned int linkIdx);
		MotorJointPtr getMotorJointPtr(const unsigned int motorJointIdx);
		const LinkPtr& getLinkPtr(const unsigned int linkIdx) const;
		const MotorJointPtr& getMotorJointPtr(const unsigned int motorJointIdx) const;
		bool isComplete() const;
		const unsigned int getNumOfLink() const;
		const unsigned int getNumOfJoint() const;

		/*!
		* \brief Link addition & deletion function
		* 
		*/
		void addLink(const LinkPtr& linkPtr);
		void deleteLink(unsigned int linkIdx);
		void insertLink(unsigned int linkIdx, const LinkPtr& linkPtr);

		/*!
		* \brief Joint addition & deletion function
		* 
		*/
		void addMotorJoint(const MotorJointPtr& motorJointPtr);
		void deleteMotorJoint(unsigned int motorJointIdx);
		void insertMotorJoint(unsigned int motorJointIdx, const MotorJointPtr& motorJointPtr);

		/*!
		* \brief Link-Joint mate function
		*
		*/
		void addMate(const unsigned int motorJointIdx, const SE3& parentT, const SE3& childT);

		/*!
		* \brief Assembly complete function
		*        this function make baseframe and
		*        calculate initial links global frame, inertia and joint screw  
		*/
		void completeAssembling();

		/*!
		* \brief make state function
		*        this function returns initial link, joint state
		*/
		StatePtr makeState() const;


	public:
		/*!
		* \brief Serial Open Chain kinematics functions
		*/
		 void solveForwardKinematics(State& state);
		 void solveDiffForwardKinematics(State& state);
		 void solve2ndDiffForwardKinematics(State& state);

		 void solveInverseKinematics(State& state, const SE3& goalT);

		 void solveJacobian(State& state);
		 void solveJacobianDot(State & state);

		 ///< calculate joint exponential(exp([S_i] * theta_i) if this value is not up to date
		 void updateJointStateExponetial(State & state, const unsigned int jointIndex);
		 void updateAccumulatedT(State & state);

	public:
		/*!
		* \brief Serial Open Chain dynamics functions
		*/
		/*!
		* \brief Inverse Dynamics function
		*        endeffectorF is represented by body frame
		*/
		void solveInverseDynamics(State& state, const dse3& endeffectorF = (dse3::Zero()));
		void solveFowardDynamics(State& state);

		MatrixX differentiateInverseDynamics(State& state, const MatrixX& dqdp, const MatrixX& dqdotdp, const MatrixX& dqddotdp);

	};


	class Mate
	{
	private:
		/*!
		* \brief Mate class member variable
		*/
		MotorJointPtr _motorJoint;
		// parent link(prior link) SE3 w.r.t joint frame
		SE3 _parentT;
		// child link(posterior link) SE3 w.r.t joint frame
		SE3 _childT;

	public:
		/*!
		* \brief Mate class member functions
		*/
		// constructor & destructor
		Mate(const MotorJointPtr& motorJoint = MotorJointPtr(),
			const SE3& parentT = SE3(), const SE3& childT = SE3());
		~Mate();

		// set-functions
		void setMotorJoint(const MotorJointPtr& motorJoint);
		void setParentT(const SE3& parentT);
		void setChildT(const SE3& childT);

		// get-functions
		const MotorJointPtr& getMotorJoint() const;
		const SE3& getParentT() const;
		const SE3& getChildT() const;
	};


	class SerialOpenChainLink
	{
	private:
		SE3 _M; ///< initial link SE3 w.r.t base(fixed) frame
		Inertia _G; ///< initial link inertia w.r.t base(fixed) frame 
	public:
		/*!
		* \brief SerialOpenChainLink class member functions
		*/
		// constructor & destructor
		SerialOpenChainLink(const SE3& M = SE3(), const Inertia& G = Inertia());
		~SerialOpenChainLink();

		// set-functions
		void setM(const SE3& M);
		void setG(const Inertia& G);

		// get-functions
		const SE3& getM() const;
		const Inertia& getG() const;
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	};


	class SerialOpenChainJoint
	{
	private:
		se3 _screw; ///< initial joint screw w.r.t base(fixed) frame
	public:
		/*!
		* \brief SerialOpenChainMate class member functions
		*/
		// constructor & destructor
		SerialOpenChainJoint();
		SerialOpenChainJoint(const se3& screw);
		~SerialOpenChainJoint();

		// set-functions
		void setScrew(const se3& screw);

		// get-functions
		const se3& getScrew() const;
		
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	};
}