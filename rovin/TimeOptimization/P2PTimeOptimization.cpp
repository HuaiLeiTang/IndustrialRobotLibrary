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
		std::vector<unsigned int> allMVCPointsFlag;
		_tree->_avp_rrt->_topp->calculateAllMVCPoint();
		allMVCPoints = _tree->_avp_rrt->_topp->getAllMVCPoint();
		allMVCPointsFlag = _tree->_avp_rrt->_topp->getAllMVCPointFlag();
		LOG("calculate allMVCPoints complete.");

		std::vector<SwitchPoint> allSwitchPoint;
		_tree->_avp_rrt->_topp->calculateAllSwitchPoint();
		allSwitchPoint = _tree->_avp_rrt->_topp->getAllSwitchPoint();
		LOG("calculate allSwitchPoint complete.");

		// calcuate all LC ----> 이 알고리즘이 너무 비효율적인거 같은데.. 다른 아이디어?
		bool A1switch;
		std::vector<std::list<Vector2, Eigen::aligned_allocator<Vector2>>> LC;
		A1switch = calculateLimitingCurves(allMVCPoints, allMVCPointsFlag, allSwitchPoint, LC);
		LOG("calculate Limiting curves complete.");

		// save LC
		for (unsigned int i = 0; i < LC.size(); i++)
		{
			std::string s_string = "C:/Users/crazy/Desktop/Time optimization/LC/s_LC";
			std::string sdot_string = "C:/Users/crazy/Desktop/Time optimization/LC/sdot_LC";
			s_string = s_string + std::to_string(i) + ".txt";
			sdot_string = sdot_string + std::to_string(i) + ".txt";
			saveLC(LC[i], s_string, sdot_string);
		}

		// case A1
		if (A1switch == false)
		{
			LOG("A1 case, Limiting curve reaches sdot = 0. AVP failure.");
			return false;
		}

		// calculate CLC
		std::vector<Vector2, Eigen::aligned_allocator<Vector2>> CLC;
		calculateCLC(LC, CLC);

		// save MVC
		saveMVC(allMVCPoints);

		// save CLC
		std::vector<Real> s_CLC;
		std::vector<Real> sdot_CLC;
		for (unsigned int i = 0; i < CLC.size(); i++)
		{
			s_CLC.push_back(CLC[i](0));
			sdot_CLC.push_back(CLC[i](1));
		}
		saveRealVector2txt(s_CLC, "C:/Users/crazy/Desktop/Time optimization/s_CLC.txt");
		saveRealVector2txt(sdot_CLC, "C:/Users/crazy/Desktop/Time optimization/sdot_CLC.txt");
		
		

		// save switch point
		std::vector<Real> s_sw;
		std::vector<Real> sdot_sw;
		for (unsigned int i = 0; i < allSwitchPoint.size(); i++)
		{
			//std::cout << allSwitchPoint[i]._id << std::endl;
			s_sw.push_back(allSwitchPoint[i]._s);
			sdot_sw.push_back(allSwitchPoint[i]._sdot);
		}
		saveRealVector2txt(s_sw, "C:/Users/crazy/Desktop/Time optimization/s_sw.txt");
		saveRealVector2txt(sdot_sw, "C:/Users/crazy/Desktop/Time optimization/sdot_sw.txt");

		return true;
	}

	bool Extend::calculateLimitingCurves(const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints, 
		const std::vector<unsigned int>& allMVCPointsFlag, const std::vector<SwitchPoint>& allSwitchPoint,
		std::vector<std::list<Vector2, Eigen::aligned_allocator<Vector2>>>& LC)
	{
		Real s_cur, sdot_cur, s_swi, sdot_swi, sdot_MVC, s_next, sdot_next, slope;
		Real ds = _tree->_avp_rrt->_topp->getStepSize();
		Real si = _tree->_avp_rrt->_topp->getInitialParam();
		Real sf = _tree->_avp_rrt->_topp->getFinalParam();
		int flag;
		unsigned int idx;
		unsigned int numOfswi = allSwitchPoint.size();


		Vector2 alphabeta;
		Real alpha_cur, beta_cur;
		std::list<Vector2, Eigen::aligned_allocator<Vector2>> tmp_list;

		for (unsigned int i = 0; i < numOfswi; i++)
		{
			tmp_list.clear();
			LC.push_back(tmp_list);

			s_swi = round(allSwitchPoint[i]._s / ds) * ds;
			sdot_swi = allSwitchPoint[i]._sdot;

			// backward integration
			s_cur = s_swi;
			sdot_cur = sdot_swi;
			LC[i].push_front(Vector2(s_cur, sdot_cur));

			if (allSwitchPoint[i]._id == SwitchPoint::SINGULAR)
			{
				Real lambda = allSwitchPoint[i]._lambda;
				unsigned int numOfSPInt = 3;
				for (unsigned int j = 0; j < numOfSPInt; j++)
				{
					s_cur -= ds;
					sdot_cur -= lambda * ds;
					LC[i].push_front(Vector2(s_cur,sdot_cur));
				}
			}

			while (s_cur > si)
			{
				alpha_cur = _tree->_avp_rrt->_topp->determineAlphaBeta(s_cur, sdot_cur)(0);
				_tree->_avp_rrt->_topp->backwardIntegrate(s_cur, sdot_cur, alpha_cur);
				
				idx = round(s_cur / ds);
				sdot_MVC = allMVCPoints[idx](1);
				//sdot_MVC2 = _tree->_avp_rrt->_topp->calculateMVCPoint(s_cur, flag);

				flag = allMVCPointsFlag[idx];

				if (sdot_cur < 1e-4)
					return false;

				LC[i].push_front(Vector2(s_cur, sdot_cur));
				if (sdot_cur > sdot_MVC)
				{
					//sdot_cur = sdot_MVC;
					//LC[i].pop_front(); LC[i].push_front(Vector2(s_cur, sdot_cur));
					//break;
					if (flag == 2)
					{
						//LC[i].pop_front(); LC[i].push_front(Vector2(s_cur, sdot_cur));
						s_cur += ds;
						LC[i].pop_front();
						break;
					}
					else if (flag == 1)
					{
						s_cur += ds;
						LC[i].pop_front();
						while (true)
						{
							s_next = s_cur - ds;
							if (s_next < si)
								break;

							idx = round(s_next / ds);
							sdot_next = allMVCPoints[idx](1);
							flag = allMVCPointsFlag[idx];

							// 가속도/토크 constraint 만나면 나가라!!
							if (flag == 2)
								break;

							alphabeta = _tree->_avp_rrt->_topp->determineAlphaBeta(s_cur, sdot_cur);
							alpha_cur = alphabeta(0);
							beta_cur = alphabeta(1);

							slope = (sdot_next - sdot_cur) / (s_next - s_cur);

							s_cur = s_next; sdot_cur = sdot_next;
							LC[i].push_back(Vector2(s_cur, sdot_cur));

							//if (beta_cur >= slope && slope >= alpha_cur)
							//{
							//	s_cur = s_next; sdot_cur = sdot_next;
							//	
							//}



						}
					}
					
				}
			}

			// forward integration
			s_cur = s_swi;
			sdot_cur = sdot_swi;

			if (allSwitchPoint[i]._id == SwitchPoint::SINGULAR)
			{
				Real lambda = allSwitchPoint[i]._lambda;
				unsigned int numOfSPInt = 3;
				for (unsigned int j = 0; j < numOfSPInt; j++)
				{
					s_cur += ds;
					sdot_cur += lambda * ds;
					LC[i].push_back(Vector2(s_cur,sdot_cur));
				}
			}

			while (s_cur < sf)
			{
				beta_cur = _tree->_avp_rrt->_topp->determineAlphaBeta(s_cur, sdot_cur)(1);
				_tree->_avp_rrt->_topp->forwardIntegrate(s_cur, sdot_cur, beta_cur);
				//sdot_MVC = _tree->_avp_rrt->_topp->calculateMVCPoint(s_cur, flag);
				idx = round(s_cur / ds);
				sdot_MVC = allMVCPoints[idx](1);
				flag = allMVCPointsFlag[idx];

				if (sdot_cur < 1e-4)
					return false;
				LC[i].push_back(Vector2(s_cur, sdot_cur));
				if (sdot_cur > sdot_MVC)
				{
					//sdot_cur = sdot_MVC;
					if (flag == 2)
					{
						//LC[i].pop_back(); LC[i].push_back(Vector2(s_cur, sdot_cur));
						s_cur -= ds;
						LC[i].pop_back();
						break;
					}
					else if (flag == 1)
					{
						s_cur -= ds;
						LC[i].pop_back();
						while (true)
						{
							s_next = s_cur + ds;
							idx = round(s_next / ds);
							sdot_next = allMVCPoints[idx](1);
							flag = allMVCPointsFlag[idx];

							// 가속도/토크 constraint 만나면 나가라!!
							if (flag == 2)
								break;

							alphabeta = _tree->_avp_rrt->_topp->determineAlphaBeta(s_cur, sdot_cur);
							alpha_cur = alphabeta(0);
							beta_cur = alphabeta(1);

							slope = (sdot_next - sdot_cur) / (s_next - s_cur);

							//
							s_cur = s_next; sdot_cur = sdot_next;
							LC[i].push_back(Vector2(s_cur, sdot_cur));

							//if (beta_cur >= slope && slope >= alpha_cur)
							//{
							//	s_cur = s_next; sdot_cur = sdot_next;
							//	LC[i].push_back(Vector2(s_cur, sdot_cur));
							//}
							//else if (beta_cur < slope)
							//{
							//	// forward integrate..,.
							//}



						}
					}
				}
			}
		}
		return true;
	}

	void Extend::calculateCLC(std::vector<std::list<Vector2, Eigen::aligned_allocator<Vector2>>>& LC,
		std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& CLC)
	{
		unsigned int LCsize = LC.size();

		std::vector<std::list<Vector2, Eigen::aligned_allocator<Vector2>>::iterator> it_begin;
		std::vector<std::list<Vector2, Eigen::aligned_allocator<Vector2>>::iterator> it_end;
		
		Real s_begin, s_end;
		
		for (unsigned int i = 0; i < LCsize; i++)
		{
			std::list<Vector2, Eigen::aligned_allocator<Vector2>>::iterator it_tmp;
			it_tmp = LC[i].end(); it_tmp--;
			it_end.push_back(it_tmp);
		}

		// find s_end
		Real tmp = -std::numeric_limits<Real>::max();
		for (unsigned int i = 0; i < LCsize; i++)
		{
			if ((*(it_end[i]))(0) > tmp)
				tmp = (*(it_end[i]))(0);
		}
		s_end = tmp;

		for (unsigned int i = 0; i < LCsize; i++)
		{
			std::list<Vector2, Eigen::aligned_allocator<Vector2>>::iterator it_tmp;
			it_tmp = LC[i].begin();
			it_begin.push_back(it_tmp);
		}

		// find s_start
		tmp = std::numeric_limits<Real>::max();
		for (unsigned int i = 0; i < LCsize; i++)
		{
			if ((*(it_begin[i]))(0) < tmp)
				tmp = (*(it_begin[i]))(0);
		}
		s_begin = tmp;

		// find CLC
		Real s_cur = s_begin;
		std::vector<unsigned int> idx;
		Real ds = _tree->_avp_rrt->_topp->getStepSize();
		while (s_cur < s_end)
		{
			idx.clear();

			for (unsigned int i = 0; i < LCsize; i++)
			{
				if ((*(it_end[i]))(0) < s_cur)
				{
					it_end.erase(it_end.begin() + i);
					it_begin.erase(it_begin.begin() + i);
					LC.erase(LC.begin() + i);
					LCsize -= 1;
				}
			}

			for (unsigned int i = 0; i < LCsize; i++)
			{
				if ((*(it_begin[i]))(0) <= s_cur)
					idx.push_back(i); // 우선 비교 idx만 저장
			}

			Real tmp_min = std::numeric_limits<Real>::max();
			if (idx.size() == 0)
				CLC.push_back(Vector2(s_cur, tmp_min));
			else
			{
				
				for (unsigned int i = 0; i < idx.size(); i++)
				{
					if ((*(it_begin[idx[i]]))(1) < tmp_min)
						tmp_min = (*(it_begin[idx[i]]))(1);
					it_begin[idx[i]]++;
				}
				CLC.push_back(Vector2(s_cur, tmp_min));
			}
			s_cur += ds;
		}

		for (unsigned int i = 0; i < LCsize; i++)
		{
			if ((*(it_end[i]))(0) < s_cur)
			{
				it_end.erase(it_end.begin() + i);
				it_begin.erase(it_begin.begin() + i);
				LC.erase(LC.begin() + i);
				LCsize -= 1;
			}
		}

		for (unsigned int i = 0; i < LCsize; i++)
		{
			if ((*(it_begin[i]))(0) < s_cur)
				idx.push_back(i);
		}

		Real tmp_min = std::numeric_limits<Real>::max();
		if (idx.size() == 0)
			CLC.push_back(Vector2(s_cur, tmp_min));
		else
		{
			for (unsigned int i = 0; i < idx.size(); i++)
			{
				if ((*(it_begin[idx[i]]))(1) < tmp_min)
					tmp_min = (*(it_begin[idx[i]]))(1);
				it_begin[idx[i]]++;
			}
			CLC.push_back(Vector2(s_cur, tmp_min));
		}

		
	}
}