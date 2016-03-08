#include "P2PTimeOptimization.h"

namespace rovin
{
	/////////////////// AVP_RRT class ///////////////////

	AVP_RRT::AVP_RRT() : _robot(SerialOpenChainPtr()) {}

	AVP_RRT::~AVP_RRT()
	{
		unsigned int treesize = _tree.size();
		for (unsigned int i = 0; i < treesize; i++)
			delete _tree[i];
		delete _topp;
	}
	
	AVP_RRT::AVP_RRT(const SerialOpenChainPtr & robot, CONSTRAINT_TYPE constraintType)
	{
		_robot = robot;
		_dof = robot->getNumOfJoint();

		Real ds = 1e-3, vi = 0, vf = 0, si = 0, sf = 1;
		_topp = new TOPP(_robot, vi, vf, ds, si, sf, constraintType);
	}

	void AVP_RRT::initialization()
	{

	}

	void AVP_RRT::makeRandomConfig(VectorX & qrand)
	{

	}

	void AVP_RRT::generateTrajectory()
	{
		initialization();

		// make tree
		unsigned int waypointsize = _waypoints.size();
		unsigned int treesize = _tree.size();
		LOGIF((treesize + 1 == waypointsize), "AVP_RRT::generateTrajectory error : waypoint size must be one more than tree size.");
		for (unsigned int i = 0; i < treesize; i++)
			_tree[i] = new Tree(this, _waypoints[i], _waypoints[i + 1]);

		//

	}


	/////////////////// Tree class ///////////////////

	Tree::Tree(AVP_RRT * avp_rrt, const WayPoint & initPoint, const WayPoint & finalPoint) : _avp_rrt(avp_rrt)
	{
		_initPoint = initPoint;
		_finalPoint = finalPoint;
		
		Vertex Vinit;
		Vinit.setconfig(_initPoint.getJointq());
		Vinit.setinterval(Vector2(0, 0));
		Vinit.setparentVertexIdx(-1);

		_vertices.push_back(Vinit);

		_extend = new Extend(this);
	}

	bool Tree::extend(VectorX & qrand)
	{
		//_extend = new Extend(this, qrand);
		_extend->nearestNeighbor(qrand);
		_extend->interpolate();
		if (_extend->collisioncheck())
		{
			if (_extend->AVP())
			{
				_extend->update();
				return true;
			}
			else
				return false;
		}
		else
			return false;

		// collisioncheck 해서 collision 되면 어디 다시 실행? --> path를 다시 만들어야 되나요??
		// AVP 에서 timeparameteriable 하지 않으면 어디 다시 실행?
	}

	void Tree::makepath()
	{
	}

	/////////////////// Extend class ///////////////////

	void Extend::update()
	{
		_Vnew.setconfig(_qnew);
		//_Vnew.setinpath(_pathnew);
		_Vnew.setinterval(_intervalnew);
		_Vnew.setparentVertexIdx(_nearestVertexIdx);
		_tree->_vertices.push_back(_Vnew);
	}

	void Extend::nearestNeighbor(VectorX& qrand)
	{
		_qrand = qrand;
		// find Vnear
		// TODO

	}

	void Extend::interpolate()
	{
		// make Pnew and qnew
		// TODO
	}

	bool Extend::collisioncheck()
	{
		return true;
	}

	bool Extend::AVP()
	{
		// calculate Vnew interval sdot_min, sdot_max
		
		// make topp
		//_tree->_avp_rrt->_topp->setJointTrajectory(_pathnew);

		// Step A. Computing the limiting curves
		// Acalculate MVC
		std::vector<Vector2, Eigen::aligned_allocator<Vector2>> allMVCPoints;
		_tree->_avp_rrt->_topp->calculateAllMVCPoint();
		allMVCPoints = _tree->_avp_rrt->_topp->getAllMVCPoint();
		LOG("calculate allMVCPoints complete.")

		// calculate CLC  ----> 이 알고리즘이 너무 비효율적인거 같은데.. 다른 아이디어?
		std::vector<SwitchPoint> allSwitchPoint;
		_tree->_avp_rrt->_topp->calculateAllSwitchPoint();
		allSwitchPoint = _tree->_avp_rrt->_topp->getAllSwitchPoint();
		LOG("calculate allSwitchPoint complete.")

		// calcuate all LC ----> 이 알고리즘이 너무 비효율적인거 같은데.. 다른 아이디어?
		// A1
		bool A1switch;
		std::vector<std::list<Real>> LC;
		A1switch = calculateLimitingCurves(allMVCPoints, allSwitchPoint, LC);
		LOG("calculate Limiting curves complete.")

		if (A1switch == false)
		{
			LOG("A1 case, Limiting curve reaches sdot = 0. AVP failure.");
			return false;
		}




		return true;
	}

	bool Extend::calculateLimitingCurves(const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints, const std::vector<SwitchPoint>& allSwitchPoint, std::vector<std::list<Real>>& LC)
	{
		Real s_cur, sdot_cur, s_swi, sdot_swi, sdot_MVC;
		Real ds = _tree->_avp_rrt->_topp->getStepSize();
		Real si = _tree->_avp_rrt->_topp->getInitialParam();
		Real sf = _tree->_avp_rrt->_topp->getFinalParam();
		int flag, idx;
		unsigned int numOfswi = allSwitchPoint.size();

		Real alpha_cur, beta_cur;

		for (unsigned int i = 0; i < numOfswi; i++)
		{
			std::list<Real> tmp_list;
			LC.push_back(tmp_list);

			// switchpoint rounding
			s_swi = round(allSwitchPoint[i]._s / ds) * ds;
			sdot_swi = allSwitchPoint[i]._sdot;

			// backward integration
			s_cur = s_swi;
			sdot_cur = sdot_swi;
			LC[i].push_front(sdot_cur);

			if (allSwitchPoint[i]._id == SwitchPoint::SINGULAR)
			{
				Real lambda = allSwitchPoint[i]._lambda;
				unsigned int numOfSPInt = 3;
				for (unsigned int i = 0; i < numOfSPInt; i++)
				{
					s_cur -= ds;
					sdot_cur -= lambda * ds;
					LC[i].push_front(sdot_cur);
				}
			}

			while (s_cur > si)
			{
				alpha_cur = _tree->_avp_rrt->_topp->determineAlphaBeta(s_cur, sdot_cur)(0);
				_tree->_avp_rrt->_topp->backwardIntegrate(s_cur, sdot_cur, alpha_cur);
				//sdot_MVC = _tree->_avp_rrt->topp->calculateMVCPoint(s_cur, flag);
				sdot_MVC = allMVCPoints[s_cur / ds](1);

				if (sdot_cur < 1e-4)
					return false;

				LC[i].push_front(sdot_cur);

				if (sdot_cur > sdot_MVC)
				{
					LC[i].pop_front(); LC[i].push_front(sdot_MVC);
					break;
				}
			}
			while (s_cur > si)
			{
				s_cur -= ds;
				sdot_cur = std::numeric_limits<Real>::max();
				LC[i].push_front(sdot_cur);
			}

			// forward integration
			s_cur = s_swi;
			sdot_cur = sdot_swi;

			if (allSwitchPoint[i]._id == SwitchPoint::SINGULAR)
			{
				Real lambda = allSwitchPoint[i]._lambda;
				unsigned int numOfSPInt = 3;
				for (unsigned int i = 0; i < numOfSPInt; i++)
				{
					s_cur += ds;
					sdot_cur += lambda * ds;
					LC[i].push_front(sdot_cur);
				}
			}

			while (s_cur > si)
			{
				beta_cur = _tree->_avp_rrt->_topp->determineAlphaBeta(s_cur, sdot_cur)(1);
				_tree->_avp_rrt->_topp->forwardIntegrate(s_cur, sdot_cur, beta_cur);
				//sdot_MVC = _tree->_avp_rrt->topp->calculateMVCPoint(s_cur, flag);
				sdot_MVC = allMVCPoints[s_cur / ds](1);

				if (sdot_cur < 1e-4)
					return false;
				LC[i].push_back(sdot_cur);

				if (sdot_cur > sdot_MVC)
				{
					LC[i].pop_back(); LC[i].push_back(sdot_MVC);
					break;
				}
			}

			while (s_cur > si)
			{
				s_cur -= ds;
				sdot_cur = std::numeric_limits<Real>::max();
				LC[i].push_back(sdot_cur);
			}
		}

		return true;
	}

	

}