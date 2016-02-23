#pragma once

#include <rovin\Dynamics\SerialOpenChain.h>
#include <rovin\Math\Interpolation.h>
#include <rovin\Optimizer\NonlinearOptimization.h>

#include <list>

namespace rovin {

	class TOPP;
	class SwitchPoint;

	class TOPP
	{


	private:
		MatrixX _q_data;
		
		Spline _q;
		Spline _dqds;
		Spline _ddqdds;

		Real _vi;
		Real _vf;

		std::list<Real> _s;
		std::list<Real> _sdot;

		VectorX _torqueConstraint;
		
		Real _ds;

		
		std::vector<SwitchPoint> _switchPoint;
		
		//std::vector<Vector2> _sdotTrajectory;

		SerialOpenChainPtr _soc;

		unsigned int _dof;

		NonlinearOptimization _nop;

	public:
		TOPP(const MatrixX& q_data, const SerialOpenChainPtr& soc, const Real ds, const Real vi, const Real vf);

		bool checkMVCCondition(Real alpha, Real beta);
		
		VectorX& calculateA(Real s);
		std::vector<VectorX>& calculateBandC(Real s, Real sdot);

		bool checkDynamicSingularity(const VectorX& a);
		Vector2& determineAlphaBeta(Real s, Real sdot);

		void farwardIntegrate(Real& s, Real& sdot, Real sddot);
		void backwardIntegrate(Real& s, Real& sdot, Real sddot);
		
		void findNearestSwitchPoint(Real s, Real sdot);
		Real calulateMVCPoint(Real s);

		void generateTrajectory();
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