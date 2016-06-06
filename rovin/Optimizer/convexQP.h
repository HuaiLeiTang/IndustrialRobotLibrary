#pragma once

#include "rovin\Math\Common.h"
#include "rovin\Math\Constant.h"

using namespace rovin;

namespace irLib
{
	namespace irOpt
	{
		class convexQP;


		class convexQP
		{
		public:

		public:
			static double solveEqualityConstraintedQP(const MatrixX& G, const VectorX& c, const MatrixX& A, const VectorX& b);

		
		};
	}
}
