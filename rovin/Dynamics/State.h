/*!
 *	\file	State.h
 *	\date	2016.01.22
 *	\author	Youngsuk (crazyhys@gmail.com)
 *	\brief	State class
 *          this class has link state & joint state
*/

#pragma once

#include <memory>
#include <vector>

#include <rovin\Math\LieGroup.h>


namespace rovin {

	class State;
	class JointState;
	class LinkState;
	
	typedef std::shared_ptr<State> StatePtr;
	typedef std::shared_ptr<JointState> JointStatePtr;
	typedef std::shared_ptr<LinkState> LinkStatePtr;

	//////////////////////////////////////////////////////////////
	//						 STATE CLASS			     		//
	//////////////////////////////////////////////////////////////
	class State
	{	
		/*!
		* \brief Information about whether each state is up to date or not
		*/
		const enum STATE_INFO
		{
			LINKS_POS = 1 << 0,		//	Position(SE3) of all links is up to date
			LINKS_VEL = 1 << 1,		//	Velocity(se3) of all links is up to date
			LINKS_ACC = 1 << 2,		//	Accelaration(se3) of all links is up to date
			JOINTS_T_FROM_BASE = 1 << 3,	//	Product of transform from base to each joint is up to date
			JOINTS_JACOBIAN = 1 << 4,		//	Jacobian(theta) is up to date
			JOINTS_JACOBIAN_DOT = 1 << 5,	//	(dJacobian(theta) / dt) is up to date

			ALL_LINKS = LINKS_POS | LINKS_VEL | LINKS_ACC,
			ALL_JOINTS = JOINTS_T_FROM_BASE | JOINTS_JACOBIAN | JOINTS_JACOBIAN_DOT,
			ALL_INFO = LINKS_POS | LINKS_VEL | LINKS_ACC | JOINTS_T_FROM_BASE | JOINTS_JACOBIAN | JOINTS_JACOBIAN_DOT,
		};

	private:
		/*!
		* \brief State class member variables
		*/
		// State total degree of freedom
		unsigned int _dof;

		// 'LinkState's are stored in assembeld order
		std::vector< LinkState, Eigen::aligned_allocator< LinkState >> _linkState;

		//	'JointState's are stored in assembled order.
		std::vector< JointState, Eigen::aligned_allocator< JointState >> _jointState;
		
		//	i-th bit of '_stateInfoUpToDate' contains up-to-date state of i-th 'STATE_INFO'
		int	_stateInfoUpToDate;

	public:
		/*!
		* \brief State class member functions
		*/

		// constructor & desctructor
		State(unsigned int dof);
		~State();

		// get-functions
		const int getDof() const;
		const JointState& getJointState(const unsigned int jointIndex);
		const LinkState& getLintState(const unsigned int linkIndex);

		// set-functions
		void setJointStatePos(JointState& jointstate, const Real q);
		void setJointStatePos(const unsigned int jointIndex, const Real q);
		void setJointStateVel(JointState& jointstate, const Real qdot);
		void setJointStateVel(const unsigned int jointIndex, const Real qdot);
		void setJointStateAcc(JointState& jointstate, const Real qddot);
		void setJointStateAcc(const unsigned int jointIndex, const Real qddot);

		// joint state(pos, vel, acc) add functions
		void addJointStatePos(JointState& jointstate, const Real q);
		void addJointStatePos(const unsigned int jointIndex, const Real q);
		void addJointStateVel(JointState& jointstate, const Real qdot);
		void addJointStateVel(const unsigned int jointIndex, const Real qdot);
		void addJointStateAcc(JointState& jointstate, const Real qddot);
		void addJointStateAcc(const unsigned int jointIndex, const Real qddot);

		// other functions
		/*!
		* checkJointInfoUpToData function description
		* 
		*/
		bool checkStateInfoUpToDate(int infoIdx);
		/*!
		* updataJointInfoUpToData funcstion description
		* 
		*/
		void updateStateInfoUpToDate(int infoIdx, bool upToDate = true);

	};

	//////////////////////////////////////////////////////////////
	//						LINK STATE CLASS					//
	//////////////////////////////////////////////////////////////
	class LinkState
	{
	private:
		/*!
		* \brief LinkState class member variables
		*/
		// Link frame SE3 w.r.t base(fixed) frame
		SE3 _T;

		// Link frame generalized velocity w.r.t base(fixed) frame
		se3 _V;

		// Link frame generalized velocity derivative w.r.t base(fixed) frame
		se3 _VDot;

		// Articulated inertia of link which is used in forward dynamics (어디 프레임에서 기술된건지 알아볼 것)
		Matrix6 _Ja;

		// Bias force which is also used in forward dynamics (어디 프레임에서 기술된건지 알아볼 것)
		dse3 _b;

	public:
		/*!
		* \brief LinkState class member variables
		*/
		// constructor & destructor
		LinkState();
		~LinkState();

		// set-function
		void setLinkSE3(const SE3& T);
		void setLinkVel(const se3& V);
		void setLinkVelDot(const se3& VDot);
		void setLinkArtInertia(const Matrix6& Ja);
		void setLinkBiasforce(const dse3& b);

		// get-function
		const SE3& getLinkSE3() const;
		const se3& getLinkVel() const;
		const se3& getLinkVelDot() const;
		const Matrix6& getLinkArtInertia() const;
		const dse3& getLinkBiasforce() const;

		// deep-copy
		LinkStatePtr copy() const;

	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	};


	//////////////////////////////////////////////////////////////
	//						JOINT STATE CLASS					//
	//////////////////////////////////////////////////////////////
	class JointState
	{
		/*!
		* \brief Information about whether each joint state is up to date or not
		*/
		const enum JOINT_INFO
		{
			TRANSFORM = 1 << 0, // exp([S_i] * theta_i)(variable _T) is up to date
			JACOBIAN = 1 << 1, // need description
			JACOBIAN_DOT = 1 << 2, // need description
			ALL_INFO = TRANSFORM | JACOBIAN | JACOBIAN_DOT,
		};

	private:
		/*!
		* \brief LinkState class member variables
		*/
		// torque
		Real _tau;

		// exp([S_i] * theta_i)
		SE3 _T;

		// Products of exp([S_i] * theta_i) from base
		SE3 _accumulatedT;

		// screw w.r.t base(fixed) frame
		se3 _screw;

		// screw derivative w.r.t base(fixed) frame
		se3 _screwDot;

		// constraint force (어느 프레임에서 기술된 건지 확인 - 아마 글로벌 프레임)
		dse3 _constraintF;

		// joint angle, velocity and acceleration
		Real _q;
		Real _qdot;
		Real _qddot;

		// i - th bit of '_jointInfoUpToDate' contains up - to - date state of i - th 'JOINT_INFO'
		int _jointInfoUpToDate;

	public:
		/*!
		* \brief LinkState class member functions
		*/
		
		// constructor & destructor
		JointState();
		~JointState();

		// set-functions
		void setJointTorque(const Real tau);
		void setJointT(const SE3& T);
		void setJointAccumulatedT(const SE3& accumulatedT);
		void setJointScrew(const se3& screw);
		void setJointScrewDot(const se3& screwdot);
		void setJointConstraintF(const dse3& constraintF);
		void setJointPos(const Real q);
		void setJointVel(const Real qdot);
		void setJointAcc(const Real qddot);

		// joint state(pos, vel, acc) add functions
		void addJointPos(const Real q);
		void addJointVel(const Real qdot);
		void addJointAcc(const Real qddot);

		// get-functions
		const Real getJointTorque() const;
		const SE3& getJointT() const;
		const SE3& getJointAccumulatedT() const;
		const se3& getJointScrew() const;
		const se3& getJointScrewDot() const;
		const dse3& getJointConstraintF() const;
		const Real getJointPos() const;
		const Real getJointVel() const;
		const Real getJointAcc() const;

		// other functions
		/*!
		* checkJointInfoUpToData function description
		* it is the same as in State class function
		*/
		bool checkJointInfoUpToDate(int infoIdx);
		/*!
		* updataJointInfoUpToData funcstion description
		* it is the same as in State class function
		*/
		void updateJointInfoUpToDate(int infoIdx, bool upToDate = true);

		// deep-copy
		JointStatePtr copy() const;

	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	};

}