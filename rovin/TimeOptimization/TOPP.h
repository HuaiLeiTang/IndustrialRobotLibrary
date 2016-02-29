#pragma once

#include <rovin\Dynamics\SerialOpenChain.h>
#include <rovin\Math\Interpolation.h>
#include <rovin\Math\GaussianQuadrature.h>

#include <list>

//// temp
#include <fstream>

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
		MatrixX _q_data;
		
		BSpline<-1, -1, -1> _q;
		BSpline<-1, -1, -1> _dqds;
		BSpline<-1, -1, -1> _ddqdds;

		Real _vi, _vf; ///< initial & final velocity

		Real _si, _sf; ///< initial & final parameters

		std::list<Real> _s;
		std::list<Real> _sdot;

		VectorX _torqueConstraint;
		VectorX _velConstraint;
		VectorX _accConstraint;
		CONSTRAINT_TYPE _constraintType;
		
		Real _ds; ///< parameter step size

		std::vector<SwitchPoint> _switchPoint;

		SerialOpenChainPtr _soc;
		unsigned int _dof;

		Vector2 _minmax; ///< velocity minimum maximum (minmax(0) = min, minmix(1) = max)
		int _integrationType;

		// TOPP Result
		Real _tf_result;
		MatrixX _torque_result;

	public:
		TOPP(const MatrixX& q_data, const SerialOpenChainPtr& soc, const Real ds, 
			const Real vi, const Real vf, const Real si, const Real sf, CONSTRAINT_TYPE constrainttype);

		bool checkMVCCondition(Real alpha, Real beta);

		VectorX calculateA(Real s);
		std::vector<VectorX> calculateBandC(Real s);

		Vector2 determineAlphaBeta(Real s, Real sdot);
		void determineVelminmax(Real s);

		void farwardIntegrate(Real& s, Real& sdot, Real sddot);
		void backwardIntegrate(Real& s, Real& sdot, Real sddot);
		
		bool findNearestSwitchPoint(Real s);
		Real calculateMVCPoint(Real s);
		Real calculateMVCPointExclude(Real s, int iExclude);

		void calculateFinalTime();
		void calculateTorqueTrajectory();

		void generateTrajectory();

		// get function
		const std::list<Real>& gets() const;
		const std::list<Real>& getsdot() const;
		const Real getFinalTime() const;
		const MatrixX& getTorqueTrajectory() const;
		const unsigned int getdof() const;

	public:
		// test for MVC curve and switching point
		void calcMVC();
		void calcSPs();
		void saveRealVector2txt(std::vector<Real> in, std::string filename);
		void saveMVCandSP2txt();

		std::vector<Real> s_MVC_jk;
		std::vector<Real> sd_MVC_jk;
		std::vector<Real> s_SW_jk;
		std::vector<Real> sd_SW_jk;
		std::vector<Real>  SPID_SW_jk;

		std::vector<Real> s_FI_jk;
		std::vector<Real> sd_FI_jk;

		std::vector<Real> s_BI_jk;
		std::vector<Real> sd_BI_jk;

		std::vector<Real> s_jk;
		std::vector<Real> sdot_jk;
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