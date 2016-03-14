#include "P2PTimeOptTemp.h"

namespace rovin
{



	void AVP_RRT::generateTrajectory()
	{
		// initialization! -> way point setting

		_numSegment = _waypoints.size() - 1;
		_segmentPath.resize(_numSegment);


		_curSegment = 0;
		bool retFlag;
		do
		{
			retFlag = generateTrajectorySegment(_curSegment++);
		} while (retFlag == SUCCESS && _curSegment<_numSegment);


		std::list<VectorX>::iterator it;
		unsigned int numCol = 0;
		for (unsigned int i = 0; i < _numSegment; i++)
			numCol += _segmentPath[i].size();

		_finalPath.resize(_dof, numCol);

		int iterCol = 0;
		if (retFlag == SUCCESS)
		{
			// CONCATENATE PATH SEGMENT and SAVE AT _finalPath
			for (unsigned int i = 0; i < _numSegment; i++)
				for (it = _segmentPath[i].begin(); it != _segmentPath[i].end(); it++)
					_finalPath.col(iterCol++) = (*it);
		}



	}

	AVP_RRT::RETURNFLAG AVP_RRT::generateTrajectorySegment(int idx)
	{
		// start tree, goal tree -> waypont 로 initialize
		treeInitialization(idx);

		VectorX qrand(_dof);
		RETURNFLAG flag;

		bool connected = false;

		Vertex * candiVertex;
		Vertex * oppVertex;
		Vertex * conVertex;

		// RUN AVP-RRT!!!
		// bi-directional RRT
		int iter = 0;
		do // MAIN LOOP
		{
			//std::cout << iter << std::endl;
			makeRandomConfig(qrand);
			//if collisionchekc(qrand)
			//	continue;


			// extend from start tree to random configuration
			candiVertex = NULL;
			candiVertex = extendTree(&_startTree, qrand, true);
			if (candiVertex != NULL)
			{
				// add vertex
				_startTree.addVertex(candiVertex);
				
				// connection test to the x
				conVertex = new Vertex();
				oppVertex = new Vertex();
				if (testConnection(candiVertex, &_goalTree, true, &conVertex, &oppVertex))
				{
					// oppVertex and conVertex have same configuration but conVertex is in startTree and oppVertex is in goalTree
										
					extractPath(conVertex, oppVertex, idx);
					flag = SUCCESS;
					break;
				}
			}


			// extend from goal tree to random configuration
			candiVertex = NULL;
			candiVertex = extendTree(&_goalTree, qrand, false);
			if (candiVertex != NULL)
			{
				// add vertex
				_goalTree.addVertex(candiVertex);

				// connection test to the _startTree
				conVertex = new Vertex();
				oppVertex = new Vertex();
				if (testConnection(candiVertex, &_startTree, false, &conVertex, &oppVertex))
				{
					extractPath(oppVertex, conVertex, idx);
					flag = SUCCESS;
					break;
				}
			}

		} while (iter++ < AVP_RRT_MAX_ITER);




		if (flag == SUCCESS)
		{
			_wayPointInterval = _goalTree._nodes[0]->_interval;
		}


		_startTree.clearTree();
		_goalTree.clearTree();
		return flag;

	}


	void AVP_RRT::treeInitialization(int idx)
	{
		Vector2 tmp(0, std::numeric_limits<Real>::max());

		Vertex * startVertex = new Vertex(_waypoints[idx].getJointq(), tmp, std::list<VectorX>(), NULL);
		Vertex * goalVertex = new Vertex(_waypoints[idx + 1].getJointq(), tmp, std::list<VectorX>(), NULL);
		
		_startTree.initializeTree(startVertex);
		_goalTree.initializeTree(goalVertex);

	}

	void AVP_RRT::makeRandomConfig(VectorX & qrand)
	{
		Real lb, ub;
		//make random config and then save at qrand;
		for (unsigned int i = 0; i < _dof; i++)
		{
			lb = _robot->getMotorJointPtr(i)->getLimitPosLower();
			ub = _robot->getMotorJointPtr(i)->getLimitPosUpper();
			qrand(i) = (Real)( ((ub - lb)*rand()) / (RAND_MAX + 0.0) + lb);
		}
	}

	Vertex * AVP_RRT::extendTree(Tree * tree, const VectorX qrand, bool atStartTree) // 반환 Vertex* 로??
	{
		// find nearest vertex
		double tmpDist;
		Vertex * nVertex = tree->findNearestNeighbor(qrand, &tmpDist);




		// interpolate between nVertex.config and qnew(configuration far away from nVertex.config about stepsize)
		std::list<VectorX> Pnew;
		VectorX qnew;
		interpolate(nVertex, qrand, tmpDist, Pnew, qnew);

		// run AVP algorithm
		Vector2 endInterval;
		bool isSucceeded;
		if (atStartTree)
			isSucceeded = runAVP(Pnew, nVertex->_interval, endInterval);
		else
			isSucceeded = runAVPbackward(Pnew, nVertex->_interval, endInterval);

		// if nearest vertex is a root vertex
		if (nVertex->_parentVertex == NULL)
		{

		}

		// feasibility check! (if fails -> delete nVertex)
		// collision check
		std::list<VectorX>::iterator Piter;
		VectorX tmpVec;
		for (Piter = Pnew.begin(); Piter != Pnew.end(); Piter++)
		{
			// TO DO:
			// collision check! of Pnew.col(i)
			tmpVec = (*Piter);
			// if coliision detected, isSucceeded = false and break!
		}

		if (isSucceeded)
			return new Vertex(qnew, endInterval, Pnew, nVertex);
		else
			return NULL;

	}

	void AVP_RRT::interpolate(Vertex * nVertex, const VectorX qrand, const double dist, std::list<VectorX> & Pnew, VectorX & qnew)
	{
		if (dist < _stepsize)
			qnew = qrand;
		else
			qnew = nVertex->_config + (qrand - nVertex->_config)*_stepsize / dist;

	
		Pnew.push_front(nVertex->_config);
		Pnew.push_back(qnew);

		// TO DO: INTERPOLATE
		// from nVertex->_config to qnew!!!

	}

	bool AVP_RRT::testConnection(Vertex * vertex, Tree * tree, bool forward, /* OUTPUT */ Vertex ** cVertex, Vertex ** oVertex)
	{
		// vertex and cVertex do not belong to the input tree (belongs to the opposite tree)
		// oVertex belongs to the input tree
		// cVertex and oVertex has the same configuration, but connected to different trees

		double tmpDist;
		(*oVertex) = tree->findNearestNeighbor(vertex->_config, &tmpDist);
		//std::cout << tmpDist << std::endl;
		// connectivity test
		if (tmpDist < _stepsize)
		{
			std::list<VectorX> Ptmp;
			VectorX qtmp;
			interpolate(vertex, (*oVertex)->_config, tmpDist, Ptmp, qtmp);

			Vector2 endInterval;
			if(forward) // vertex belongs to start tree, and try to connect to goal tree
			{
				if (runAVP(Ptmp, vertex->_interval, endInterval))
				{
					if (checkIntersectionOfTwoIntervals(endInterval, (*oVertex)->_interval))
					{
						(*cVertex)->_config = (*oVertex)->_config;
						(*cVertex)->_inpath = Ptmp;
						(*cVertex)->_interval = endInterval;
						(*cVertex)->_parentVertex = vertex;
						tree->addVertex((*cVertex));
						return true;
					}
				}
			}
			else // vertex belongs to goal tree, and try to connect to start tree
			{
				if (runAVPbackward(Ptmp, vertex->_interval, endInterval))
				{
					if (checkIntersectionOfTwoIntervals(endInterval, (*oVertex)->_interval))
					{
						(*cVertex)->_config = (*oVertex)->_config;
						(*cVertex)->_inpath = Ptmp;
						(*cVertex)->_interval = endInterval;
						(*cVertex)->_parentVertex = vertex;
						tree->addVertex((*cVertex));
						return true;
					}
				}
			}
		}

		(*cVertex) = NULL;
		return false;
	}

	bool AVP_RRT::checkIntersectionOfTwoIntervals(const Vector2 & vec1, const Vector2 & vec2)
	{
		// if intersection of two intervals is null set -> return false, otherwise -> return true
		return std::max(vec1.minCoeff(), vec2.minCoeff()) <= std::min(vec1.maxCoeff(), vec2.maxCoeff());
	}

	void AVP_RRT::extractPath(Vertex * sVertex, Vertex * gVertex, int idx)
	{
		// input: sVertex (start tree), gVertex (goal tree)
		// output: save at _segmentPath[idx]
		std::list<VectorX> segPath;
		
		// path from start tree
		Vertex * curVertex = sVertex;
		segPath = curVertex->_inpath;
		curVertex = curVertex->_parentVertex;
		while (curVertex->_parentVertex != NULL)
		{
			curVertex->_inpath.pop_back();
			segPath.insert(segPath.begin(), curVertex->_inpath.begin(), curVertex->_inpath.end());
			curVertex = curVertex->_parentVertex;
		}

		// path from goal tree
		curVertex = gVertex;
		while (curVertex->_parentVertex != NULL)
		{
			curVertex->_inpath.pop_back();
			curVertex->_inpath.reverse();
			segPath.insert(segPath.end(), curVertex->_inpath.begin(), curVertex->_inpath.end());
			curVertex = curVertex->_parentVertex;
		}
		
		_segmentPath[idx] = segPath;
	}

	bool AVP_RRT::runAVP(std::list<VectorX>& Pnew, Vector2& nearInterval, Vector2 & endInterval)
	{
		endInterval(0) = 0;
		endInterval(1) = 1E4;
		return true;

		// topp initialization


		unsigned int RowSize = Pnew.front().size();
		unsigned int ColSize = Pnew.size();
		MatrixX q_data(RowSize, ColSize);
		unsigned int cnt = 0;
		for (std::list<VectorX>::iterator it = Pnew.begin(); it != Pnew.end(); it++)
			q_data.col(cnt++) = (*it);

		_topp->setJointTrajectory(q_data);

		// Step A. Computing the limiting curves
		// calculate MVC
		std::vector<Vector2, Eigen::aligned_allocator<Vector2>> allMVCPoints;
		std::vector<unsigned int> allMVCPointsFlag;
		_topp->calculateAllMVCPoint();
		allMVCPoints = _topp->getAllMVCPoint();
		allMVCPointsFlag = _topp->getAllMVCPointFlag();
		LOG("calculate allMVCPoints complete.");

		

		// calculate switching points
		std::vector<SwitchPoint> allSwitchPoint;
		_topp->calculateAllSwitchPoint();
		allSwitchPoint = _topp->getAllSwitchPoint();
		LOG("calculate allSwitchPoint complete.");

		// calculate all LC
		bool A1switch;
		std::vector<std::list<Vector2, Eigen::aligned_allocator<Vector2>>> LC;
		A1switch = calculateLimitingCurves(allMVCPoints, allMVCPointsFlag, allSwitchPoint, LC);
		LOG("calculate Limiting curves complete.");

		std::vector<std::list<Vector2, Eigen::aligned_allocator<Vector2>>> LC_copy;
		LC_copy = LC;

		// case A1
		if (A1switch == false)
		{
			LOG("A1 case, Limiting curve reaches sdot = 0. AVP failure.");
			//return false;
		}

		// calculate CLC
		std::vector<Vector2, Eigen::aligned_allocator<Vector2>> CLC;
		calculateCLC(LC, CLC);
		LOG("calculate CLC complete.");


		////////////////////////////////////////////////////////
		//////////////////////////  save  //////////////////////
		////////////////////////////////////////////////////////

		// save MVC
		saveMVC(allMVCPoints);

		// save switch point
		std::vector<Real> s_sw;
		std::vector<Real> sdot_sw;
		for (unsigned int i = 0; i < allSwitchPoint.size(); i++)
		{
			s_sw.push_back(allSwitchPoint[i]._s);
			sdot_sw.push_back(allSwitchPoint[i]._sdot);
		}
		saveRealVector2txt(s_sw, "C:/Users/crazy/Desktop/Time optimization/s_sw.txt");
		saveRealVector2txt(sdot_sw, "C:/Users/crazy/Desktop/Time optimization/sdot_sw.txt");

		// save LC
		for (unsigned int i = 0; i < LC_copy.size(); i++)
		{
			std::string s_string = "C:/Users/crazy/Desktop/Time optimization/LC/s_LC";
			std::string sdot_string = "C:/Users/crazy/Desktop/Time optimization/LC/sdot_LC";
			s_string = s_string + std::to_string(i) + ".txt";
			sdot_string = sdot_string + std::to_string(i) + ".txt";
			saveLC(LC_copy[i], s_string, sdot_string);
		}

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

		return false;
	}

	bool AVP_RRT::calculateLimitingCurves(const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints,
		const std::vector<unsigned int>& allMVCPointsFlag, const std::vector<SwitchPoint>& allSwitchPoint,
		std::vector<std::list<Vector2, Eigen::aligned_allocator<Vector2>>>& LC)
	{
		Real s_cur, sdot_cur, s_swi, sdot_swi, sdot_MVC, s_next, sdot_next, slope;
		Real ds = _topp->getStepSize();
		Real si = _topp->getInitialParam();
		Real sf = _topp->getFinalParam();
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
					LC[i].push_front(Vector2(s_cur, sdot_cur));
				}
			}

			while (s_cur > si)
			{
				alpha_cur = _topp->determineAlphaBeta(s_cur, sdot_cur)(0);
				_topp->backwardIntegrate(s_cur, sdot_cur, alpha_cur);

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
					
					s_cur += ds;
					LC[i].pop_back();
					break;

					//if (flag == 2)
					//{
					//	//LC[i].pop_front(); LC[i].push_front(Vector2(s_cur, sdot_cur));
					//	s_cur += ds;
					//	LC[i].pop_front();
					//	break;
					//}
					//else if (flag == 1)
					//{
					//	s_cur += ds;
					//	LC[i].pop_front();
					//	while (true)
					//	{
					//		s_next = s_cur - ds;
					//		if (s_next < si)
					//			break;

					//		idx = round(s_next / ds);
					//		sdot_next = allMVCPoints[idx](1);
					//		flag = allMVCPointsFlag[idx];

					//		// 가속도/토크 constraint 만나면 나가라!!
					//		if (flag == 2)
					//			break;

					//		alphabeta = _topp->determineAlphaBeta(s_cur, sdot_cur);
					//		alpha_cur = alphabeta(0);
					//		beta_cur = alphabeta(1);

					//		slope = (sdot_next - sdot_cur) / (s_next - s_cur);

					//		s_cur = s_next; sdot_cur = sdot_next;
					//		LC[i].push_back(Vector2(s_cur, sdot_cur));

					//		//if (beta_cur >= slope && slope >= alpha_cur)
					//		//{
					//		//	s_cur = s_next; sdot_cur = sdot_next;
					//		//	
					//		//}



					//	}
					//}

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
					LC[i].push_back(Vector2(s_cur, sdot_cur));
				}
			}

			while (s_cur < sf)
			{
				beta_cur = _topp->determineAlphaBeta(s_cur, sdot_cur)(1);
				_topp->forwardIntegrate(s_cur, sdot_cur, beta_cur);
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

							alphabeta = _topp->determineAlphaBeta(s_cur, sdot_cur);
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

	void AVP_RRT::calculateCLC(std::vector<std::list<Vector2, Eigen::aligned_allocator<Vector2>>>& LC,
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

		// find s_begin
		tmp = std::numeric_limits<Real>::max();
		for (unsigned int i = 0; i < LCsize; i++)
		{
			if ((*(it_begin[i]))(0) < tmp)
				tmp = (*(it_begin[i]))(0);
		}
		s_begin = tmp;

		std::cout << "s_begin : " << s_begin << std::endl;
		std::cout << "s_end : " << s_end << std::endl;

		// find CLC
		Real s_cur = s_begin;
		std::vector<unsigned int> idx;
		Real ds = _topp->getStepSize();
		while (s_cur < s_end)
		{
			idx.clear();

			for (unsigned int i = 0; i < LCsize; i++)
			{
				if ((*(it_end[i]))(0) <= s_cur)
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
			if ((*(it_begin[i]))(0) <= s_cur)
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

	bool AVP_RRT::runAVPbackward(std::list<VectorX>& Pnew, Vector2& nearInterval, Vector2 & endInterval)
	{
		endInterval(0) = 0;
		endInterval(1) = 1E4;
		return true;
	}

	void Tree::initializeTree(Vertex * rootVertex) 
	{

		_mpnnTree = new MultiANN(rootVertex->_config.size(), MPNN_TREE_MAXIMUM_NODES);

		//_nodes.clear();
		_nodes.push_back(rootVertex);
		//_numNodes = 1;

		double * rootConfig = new double[rootVertex->_config.size()];
		for (int i = 0; i < rootVertex->_config.size(); i++)
			rootConfig[i] = rootVertex->_config(i);
		_mpnnTree->AddPoint(rootConfig, (int*)rootVertex);
		delete[] rootConfig;
	}

	Vertex * Tree::findNearestNeighbor(VectorX targetConfig, double * distance) // output distance: distance from nn vertex
	{
		int index;

		double * tConfig = new double[targetConfig.size()];
		for (int i = 0; i < targetConfig.size(); i++)
			tConfig[i] = targetConfig(i);
		_mpnnTree->NearestNeighbor(tConfig, index, distance);
		delete[] tConfig;

		return _nodes[index]; // nearest node가 맞는가 -> 디버깅하면서 확인해보기.
	}

	// addVertex 부르기 전에 vertex 멤버변수 다 값 저장한다음에 넣기

	void Tree::addVertex(Vertex * vertex)
	{
		_nodes.push_back(vertex);

		double * tConfig = new double[vertex->_config.size()];
		for (int i = 0; i < vertex->_config.size(); i++)
			tConfig[i] = vertex->_config(i);
		_mpnnTree->AddPoint(tConfig, (int*)vertex);
		delete[] tConfig;
	}

}