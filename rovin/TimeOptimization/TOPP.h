#pragma once

#include <rovin\Dynamics\SerialOpenChain.h>
#include <rovin\Math\Interpolation.h>
#include <rovin\Math\GaussianQuadrature.h>
#include <list>

namespace rovin {

	class TOPP;
	class SwitchPoint;

	const enum CONSTRAINT_TYPE
	{
		TORQUE = 1 << 0,
		VELOCITY = 1 << 1,
		ACCELERATION = 1 << 2,

		TORQUE_VEL = TORQUE | VELOCITY,
		TORQUE_ACC = TORQUE | ACCELERATION,
		TORQUE_VEL_ACC = TORQUE | VELOCITY | ACCELERATION,
	};

	class TOPP
	{
	private:
		CONSTRAINT_TYPE _constraintType;
		SerialOpenChainPtr _soc;
		unsigned int _dof;

		MatrixX _q_data;
		
		BSpline<-1, -1, -1> _q;
		BSpline<-1, -1, -1> _dqds;
		BSpline<-1, -1, -1> _ddqdds;

		Real _vi, _vf;
		Real _si, _sf;
		Real _ds;

		std::list<Real> _s;
		std::list<Real> _sdot;
		std::vector<SwitchPoint> _switchPoint;

		VectorX _torqueConstraint;
		VectorX _velConstraint;
		VectorX _accConstraint;
		
		int	_nconstraints; ///< number of inequality constraints
		int _nconstraintsWithoutVel; ///< number of inequality constraints without velocity constraints

		Real _tf_result;
		MatrixX _torque_result;
	public:
		TOPP(const MatrixX& q_data, const SerialOpenChainPtr& soc, const Real ds,
			const Real vi, const Real vf, const Real si, const Real sf, CONSTRAINT_TYPE constrainttype);

	private:
		Vector2 determineAlphaBeta(Real s, Real sdot);
		std::vector<VectorX> calculateBandC(Real s);
		void calculateA(Real s, VectorX& a);
		void determineVelminmax(Real s, Vector2& minmax);
		void forwardIntegrate(Real& s, Real& sdot, Real sddot);
		void backwardIntegrate(Real& s, Real& sdot, Real sddot);
		
		bool findNearestSwitchPoint(Real s);
		Real calculateMVCPoint(Real s, int& flag);
		Real calculateMVCPointExclude(Real s, int iExclude, int& flag);

		void calculateFinalTime();
		void calculateTorqueTrajectory();

	public:
		void generateTrajectory();
		void initialization();

	public:
		void setConstraintType(CONSTRAINT_TYPE constraintType) { _constraintType = constraintType; }
		void setSerialOpenChain(const SerialOpenChainPtr& soc) { _soc = soc; }
		void setJointTrajectory(const MatrixX& q_data) { _q_data = q_data; }

		const std::list<Real>& gets() const { return _s; }
		const std::list<Real>& getsdot() const { return _sdot; }
		const Real TOPP::getFinalTime() const { return _tf_result; }
		const MatrixX& TOPP::getTorqueTrajectory() const { return _torque_result; }
		const unsigned int TOPP::getdof() const { return _dof; }
	};

	class SwitchPoint
	{
	public:
		enum SPID
		{
			SINGULAR, TANGENT, DISCONTIUOUS
		};
	public:
		Real _s;
		Real _sdot;
		SPID _id;
		Real _lambda;

	public:
		SwitchPoint(){}
		SwitchPoint(Real s, Real sdot, SPID id, Real lambda) : _s(s), _sdot(sdot),
			_id(id), _lambda(lambda) {}
	};

}