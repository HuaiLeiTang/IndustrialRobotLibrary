#pragma once

#include <rovin\Dynamics\SerialOpenChain.h>
#include <rovin\Math\Interpolation.h>


namespace rovin {

	class TOPP;



	class TOPP
	{
	private:
		std::vector<Real> _q_data;
		
		CubicSplinePtr _q;
		CubicSplinePtr _dqds;
		CubicSplinePtr _ddqdds;

		Real _ds;
		std::vector<Real> _s;
		std::vector<Real> _sdot;

		SerialOpenChainPtr _soc;



	public:
		TOPP();



	};

}