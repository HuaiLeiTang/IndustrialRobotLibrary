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

	public:
		TOPP(const MatrixX& q_data, const SerialOpenChainPtr& soc, const Real vi, const Real vf, 
			const Real ds = 1e-3, const Real si = 0, const Real sf = 1, CONSTRAINT_TYPE constraintType = TORQUE);
		TOPP(const SerialOpenChainPtr& soc, const Real vi = 0, const Real vf = 0, 
			const Real ds = 1e-3, const Real si = 0, const Real sf = 1, CONSTRAINT_TYPE constraintType = TORQUE);
		~TOPP() {}

		void generateTrajectory();
		void initialization();
		void calculateAllMVCPoint();
		void calculateAllSwitchPoint();

		void setJointTrajectory(const MatrixX& q_data);
		void setSerialOpenChain(const SerialOpenChainPtr& soc);
		void setInitialParameter(const Real si) { _si = si; }
		void setFinalParameter(const Real sf) { _sf = sf; }
		void setInitialVel(const Real vi) { _vi = vi; }
		void setFinalVel(const Real vf) { _vf = vf; }
		void setStepSize(const Real ds) { _ds = ds; }

		void setConstraintType(CONSTRAINT_TYPE constraintType) { _constraintType = constraintType; }
		
		const std::list<Real>& gets() const { return _s; }
		const std::list<Real>& getsdot() const { return _sdot; }
		const std::vector<Vector2>& getAllMVCPoint() const { return _allMVCPoints; }
		const std::vector<SwitchPoint>& getAllSwitchPoint() const { return _switchPoint; }
		const Real TOPP::getFinalTime() const { return _tf_result; }
		const MatrixX& TOPP::getTorqueTrajectory() const { return _torque_result; }
		const unsigned int TOPP::getdof() const { return _dof; }

	private:
		Vector2 determineAlphaBeta(Real s, Real sdot);
		std::vector<VectorX> calculateBandC(Real s);
		void settingconstraint();
		void makespline();
		void calculateA(Real s, VectorX& a);
		void determineVelminmax(Real s, Vector2& minmax);
		void forwardIntegrate(Real& s, Real& sdot, Real sddot);
		void backwardIntegrate(Real& s, Real& sdot, Real sddot);

		bool findNearestSwitchPoint(Real s);
		Real calculateMVCPoint(Real s, int& flag);
		Real calculateMVCPointExclude(Real s, int iExclude, int& flag);

		void calculateFinalTime();
		void calculateTorqueTrajectory();

	private:
		CONSTRAINT_TYPE _constraintType;
		SerialOpenChainPtr _soc;
		unsigned int _dof;

		MatrixX _q_data;
		MatrixX _torque_result;
		
		BSpline<-1, -1, -1> _q;
		BSpline<-1, -1, -1> _dqds;
		BSpline<-1, -1, -1> _ddqdds;

		Real _vi, _vf;
		Real _si, _sf;
		Real _ds;
		Real _tf_result;
		
		std::list<Real> _s;
		std::list<Real> _sdot;
		std::list<Real> _sddot;
		std::vector<SwitchPoint> _switchPoint;
		std::vector<Vector2> _allMVCPoints;

		VectorX _torqueConstraint;
		VectorX _velConstraint;
		VectorX _accConstraint;
		
		unsigned int _nconstraints;
		unsigned int _nconstraintsWithoutVel;
	};

	class SwitchPoint
	{
	public:
		enum SPID
		{
			SINGULAR, TANGENT, DISCONTIUOUS
		};
	public:
		SwitchPoint() {}
		SwitchPoint(Real s, Real sdot, SPID id, Real lambda) : _s(s), _sdot(sdot),
			_id(id), _lambda(lambda) {}
	public:
		Real _s;
		Real _sdot;
		SPID _id;
		Real _lambda;
	};

}