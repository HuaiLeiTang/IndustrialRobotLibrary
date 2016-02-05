#include "GivenPathTimeOptimization.h"

namespace rovin {
	GivenPathTimeOptimization::GivenPathTimeOptimization() {}

	GivenPathTimeOptimization::GivenPathTimeOptimization(const std::vector<SE3>& givenPath)
	{
		setGivenPath(givenPath);
	}

	GivenPathTimeOptimization::~GivenPathTimeOptimization()	{}

	VectorX GivenPathTimeOptimization::getq(const unsigned int i)
	{
		return _q[i];
	}

	VectorX GivenPathTimeOptimization::getqs(const unsigned int i)
	{
		return _qs[i];
	}

	VectorX GivenPathTimeOptimization::getqss(const unsigned int i)
	{
		return _qss[i];
	}

	void GivenPathTimeOptimization::setGivenPath(const std::vector<SE3>& givenPath)
	{
		_givenPath = givenPath;

		int pointsize = givenPath.size();
		for (int i = 0; i < pointsize; i++)
			_s[i] = i / Real(pointsize - 1);
	}

	void GivenPathTimeOptimization::solveMinimumTimeOptimization()
	{

	}

	VectorX GivenPathTimeOptimization::MinAcc(const Real & s, const Real & sdot)
	{

		return VectorX();
	}

	VectorX GivenPathTimeOptimization::MaxAcc(const Real & s, const Real & sdot)
	{
		return VectorX();
	}



}