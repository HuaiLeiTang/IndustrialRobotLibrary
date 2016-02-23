#pragma once

#include <rovin\Dynamics\SerialOpenChain.h>
#include <rovin\Math\Interpolation.h>
#include <rovin\Optimizer\NonlinearOptimization.h>

#include <list>

namespace rovin {

	class TOPP;

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

		std::vector<Vector2> _switchPoint;
		std::vector<Vector2> _sdotTrajectory;

		SerialOpenChainPtr _soc;

		unsigned int _dof;

		NonlinearOptimization _nop;

	public:
		TOPP(const MatrixX& q_data, const SerialOpenChainPtr& soc, const Real ds, const Real vi, const Real vf);

		bool checkMVCCondition(Real alpha, Real beta);
		
		VectorX& calculateA(Real s, Real sdot);
		std::vector<VectorX>& calculateBandC(Real s, Real sdot);

		bool checkDynamicSingularity(const VectorX& a);
		Vector2& determineAlphaBeta(Real s, Real sdot);

		void farwardIntegrate(Real& s, Real& sdot, Real sddot);
		void backwardIntegrate(Real& s, Real& sdot, Real sddot);
		
		Vector2& findNearestSwitchPoint(Real s, Real sdot);
		Real calulateMVCPoint(Real s);

		void generateTrajectory();

		

	};

}