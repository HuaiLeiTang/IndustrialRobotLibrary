#pragma once

#include <rovin\Dynamics\SerialOpenChain.h>
#include <rovin\Math\Interpolation.h>
#include <rovin\Optimizer\NonlinearOptimization.h>
#include <rovin\Math\GaussianQuadrature.h>

#include <list>

//// temp
#include <fstream>

namespace rovin {

	class TOPP;
	class SwitchPoint;

	class TOPP
	{
	public:
		MatrixX _q_data;
		
		BSpline<-1, -1, -1> _q;
		BSpline<-1, -1, -1> _dqds;
		BSpline<-1, -1, -1> _ddqdds;

		Real _vi;
		Real _vf;

		Real _si;
		Real _sf;

		std::list<Real> _s;
		std::list<Real> _sdot;

		VectorX _torqueConstraint;
		
		Real _ds;

		std::vector<SwitchPoint> _switchPoint;
		
		//std::vector<Vector2> _sdotTrajectory;

		SerialOpenChainPtr _soc;

		unsigned int _dof;

		//NonlinearOptimization _nop;

		// TOPP Result
		Real _tf_result;
		std::vector<VectorX> _torque_result;

	public:
		TOPP(const MatrixX& q_data, const SerialOpenChainPtr& soc, const Real ds, 
			const Real vi, const Real vf, const Real si, const Real sf);

		bool checkMVCCondition(Real alpha, Real beta);
		
		VectorX calculateA(Real s);
		std::vector<VectorX> calculateBandC(Real s);

		Vector2 determineAlphaBeta(Real s, Real sdot);

		void farwardIntegrate(Real& s, Real& sdot, Real sddot);
		void backwardIntegrate(Real& s, Real& sdot, Real sddot);
		
		bool findNearestSwitchPoint(Real s);
		Real calculateMVCPoint(Real s);
		Real calculateMVCPointExclude(Real s, int iExclude);

		void calculateFinalTime();
		void calculateTorqueTrajectory();

		void generateTrajectory();


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

		std::vector<Real> s_FI_jk;
		std::vector<Real> sd_FI_jk;
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