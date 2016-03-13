#include "P2PTimeOptTemp.h"

namespace rovin
{
	void AVP_RRT::generateTrajectory()
	{
		// initialization! -> way point setting

		_numSegment = _waypoints.size() - 1;


		for (unsigned int i = 0; i < _numSegment; i++)
			generateTrajectorySegment(i);

		// CONCATENATE PATH SEGMENT
		// EXTRACT FINAL PATH
	}

	AVP_RRT::RETURNFLAG AVP_RRT::generateTrajectorySegment(int idx)
	{
		// start tree, goal tree -> waypont �� initialize
		treeInitialization(idx);

		VectorX qrand;
		RETURNFLAG flag;

		bool connected = false;

		// RUN AVP-RRT!!!
		// bi-directional RRT
		int iter = 0;
		do // MAIN LOOP
		{
			makeRandomConfig(qrand);
			//if collisionchekc(qrand)
			//	continue;


			Vertex * candiVertex = extendTree(&_startTree, qrand, true);
			if (candiVertex != NULL)
			{
				// add vertex
				_startTree.addVertex(candiVertex);
				
				// connection test to the _goalTree
				Vertex * oppVertex = NULL;
				Vertex * conVertex = NULL;
				if (testConnection(candiVertex, &_goalTree, &oppVertex, &conVertex))
				{
					// oppVertex and conVertex have same configuration but conVertex is in startTree and oppVertex is in goalTree
										
					// extract path!!!!
					extractPath(conVertex, oppVertex, idx);
					return SUCCESS;
				}
			}



			extendTree(&_goalTree, qrand, false);

			// if extend succeed -> add vertex -> try to connect

			if (connected)
			{
				flag = SUCCESS;
				// -> compute trajectory OR save to final path
				break;
			}


		} while (iter++ < AVP_RRT_MAX_ITER);




		// SAVE AT PATH SEGMENT

		if (flag == SUCCESS)
		{
			_wayPointInterval = _goalTree._nodes[0]->_interval;
		}

		return flag;
	}

	void AVP_RRT::treeInitialization(int idx)
	{

		_startTree.clearTree();
		_goalTree.clearTree();

		Vertex * startVertex;
		if (idx == 0)
			startVertex = new Vertex(_waypoints[idx].getJointq(), Vector2().setZero(), std::list<VectorX>(), NULL);
		else
			startVertex = new Vertex(_waypoints[idx].getJointq(), _wayPointInterval, std::list<VectorX>(), NULL);

		Vertex * goalVertex = new Vertex(_waypoints[idx + 1].getJointq(), Vector2().setZero(), std::list<VectorX>(), NULL);


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

	Vertex * AVP_RRT::extendTree(Tree * tree, const VectorX qrand, bool atStartTree) // ��ȯ Vertex* ��??
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

		// TO DO: INTERPOLATE


	}

	bool AVP_RRT::testConnection(Vertex * vertex, Tree * tree, /* OUTPUT */ Vertex ** oVertex, Vertex ** cVertex)
	{
		// tree ���� vertex�� extended �Ǿ�����, vertex�� tree�� ����־ connection test..
		// oVertex �� tree���� ������ �Ǹ� vertex point �ƴϸ�... NULL�� ����
		// ����Ǹ� true, �ȵǸ� false..
		// �� �Լ������� avp test �ؾ���
		double tmpDist;
		(*oVertex) = tree->findNearestNeighbor(vertex->_config, &tmpDist);

		// connectivity test
		if (tmpDist < _stepsize)
		{
			std::list<VectorX> Ptmp;
			VectorX qtmp;
			interpolate(vertex, (*oVertex)->_config, tmpDist, Ptmp, qtmp);

			Vector2 endInterval;
			// ��... �̰͵� backward ���� �ؾ��ҵ�.... �ƿ� �ϴ� forward �κи� ����..
			if (runAVP(Ptmp, vertex->_interval, endInterval))
			{
				if (checkIntersectionOfTwoIntervals(endInterval, (*oVertex)->_interval))
				{
					(*cVertex)->_config = (*oVertex)->_config;
					(*cVertex)->_inpath = Ptmp;
					(*cVertex)->_interval = endInterval;
					(*cVertex)->_parentVertex = vertex;
					return true;
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




		MatrixX path;
		_segmentPath[idx] = path;
	}

	bool AVP_RRT::runAVP(std::list<VectorX>& Pnew, Vector2& nearInterval, Vector2 & endInterval)
	{
		settingtopp(Pnew);
		
		Real sdot_beg_min = nearInterval(0);
		Real sdot_beg_max = nearInterval(1);

		/* Step A. Computing the limiting curves */
		// variables
		std::vector<Vector2, Eigen::aligned_allocator<Vector2>> allMVCPoints;
		std::vector<Vector2, Eigen::aligned_allocator<Vector2>> CLC;
		std::vector<std::list<Vector2, Eigen::aligned_allocator<Vector2>>> LC;
		std::vector<unsigned int> allMVCPointsFlag;
		std::vector<SwitchPoint> allSwitchPoint;
		bool A1switch;

		_topp->calculateAllMVCPoint();
		_topp->calculateAllSwitchPoint();
		allMVCPoints = _topp->getAllMVCPoint();
		allMVCPointsFlag = _topp->getAllMVCPointFlag();
		allSwitchPoint = _topp->getAllSwitchPoint();
		A1switch = calculateLimitingCurves(allMVCPoints, allMVCPointsFlag, allSwitchPoint, LC);

		std::vector<std::list<Vector2, Eigen::aligned_allocator<Vector2>>> LC_copy;
		LC_copy = LC; ///< for save LC

		// case A1
		if (A1switch == false)
		{
			LOG("A1 case, Limiting curve reaches sdot = 0. AVP failure.");
			return false;
		}

		calculateCLC(LC, CLC);
		unsigned int Acase;
		Real sdot_beg_star;
		Acase = determineAresult(allMVCPoints, CLC, sdot_beg_star, FORWARD);

		/* Step B : determining the maximum final velocity */
		if (sdot_beg_min > sdot_beg_star)
		{
			LOG("sdot_beg_min > sdot_beg_star, the path is not traversable. AVP failure.");
			return false;
		}

		Real sdot_beg_max_star = std::min(sdot_beg_max, sdot_beg_star);
		Vector2 _nearInterval(sdot_beg_min, sdot_beg_max_star);
		std::vector<Vector2, Eigen::aligned_allocator<Vector2>> phi;
		unsigned int Bcase;
		
		Bcase = determineAVPBresult(allMVCPoints, CLC, sdot_beg_max_star, phi);

		if (Bcase == 1)
		{
			LOG("phi hits sdot = 0, the path is not traversable. AVP failure.");
			return false;
		}
		else if (Bcase == 2)
		{
			endInterval(1) = phi.back()(1);
		}
		else if (Bcase == 3)
		{
			if (Acase == 4 || Acase == 5)
				endInterval(1) = CLC.back()(1);
			else if (Acase == 2 || Acase == 3)
			{
				if (IS_VALID(allMVCPoints, CLC, phi, _nearInterval, allMVCPoints.back()(1)))
					endInterval(1) = allMVCPoints.back()(1);
				else
				{
					LOG("Case B3b, isValid is false, the path is not traversable. AVP failure.");
					return false;
				}
			}
		}
		else if (Bcase == 4)
		{
			if (IS_VALID(allMVCPoints, CLC, phi, _nearInterval, allMVCPoints.back()(1)))
				endInterval(1) = allMVCPoints.back()(1);
			else
			{
				LOG("Case B4, isValid is false, the path is not traversable. AVP failure.");
				return false;
			}
		}

		/* Step C : determining the minimum final velocity */
		endInterval(0) = findsdotminBybinearSearch(allMVCPoints, CLC, phi, _nearInterval, endInterval(1), FORWARD);
		std::cout << endInterval << std::endl;

		savedata(allMVCPoints, allSwitchPoint, LC_copy, CLC, phi);

		return true;
	}

	bool AVP_RRT::runAVPbackward(std::list<VectorX>& Pnew, Vector2& nearInterval, Vector2 & endInterval)
	{
		settingtopp(Pnew);

		Real sdot_end_min = nearInterval(0);
		Real sdot_end_max = nearInterval(1);

		/* Step A. Computing the limiting curves */
		// variables
		std::vector<Vector2, Eigen::aligned_allocator<Vector2>> allMVCPoints;
		std::vector<Vector2, Eigen::aligned_allocator<Vector2>> CLC;
		std::vector<std::list<Vector2, Eigen::aligned_allocator<Vector2>>> LC;
		std::vector<unsigned int> allMVCPointsFlag;
		std::vector<SwitchPoint> allSwitchPoint;
		bool A1switch;

		_topp->calculateAllMVCPoint();
		_topp->calculateAllSwitchPoint();
		allMVCPoints = _topp->getAllMVCPoint();
		allMVCPointsFlag = _topp->getAllMVCPointFlag();
		allSwitchPoint = _topp->getAllSwitchPoint();
		A1switch = calculateLimitingCurves(allMVCPoints, allMVCPointsFlag, allSwitchPoint, LC);

		std::vector<std::list<Vector2, Eigen::aligned_allocator<Vector2>>> LC_copy;
		LC_copy = LC; ///< for save LC

					  // case A1
		if (A1switch == false)
		{
			LOG("A1 case, Limiting curve reaches sdot = 0. AVP failure.");
			return false;
		}

		calculateCLC(LC, CLC);
		unsigned int Acase;
		Real sdot_end_star;
		Acase = determineAresult(allMVCPoints, CLC, sdot_end_star, BACKWARD);

		/* Step B : determining the maximum final velocity */
		if (sdot_end_min > sdot_end_star)
		{
			LOG("sdot_end_min > sdot_end_star, the path is not traversable. AVP backward failure.");
			return false;
		}

		Real sdot_end_max_star = std::min(sdot_end_max, sdot_end_star);  ///< �̰� �´°��� �𸣰���..
		Vector2 _nearInterval(sdot_end_min, sdot_end_max_star); ///< �̰� �´°��� �𸣰���..
		std::vector<Vector2, Eigen::aligned_allocator<Vector2>> phi;
		unsigned int Bcase;

		Bcase = determineAVPBackwardBresult(allMVCPoints, CLC, sdot_end_max_star, phi);
		std::cout << "Bcase : " << Bcase << std::endl;

		if (Bcase == 1)
		{
			LOG("phi hits sdot = 0, the path is not traversable. AVP backward failure.");
			return false;
		}
		else if (Bcase == 2)
		{
			endInterval(1) = phi.front()(1);
		}
		else if (Bcase == 3)
		{
			if (Acase == 4 || Acase == 5)
				endInterval(1) = CLC.front()(1);
			else if (Acase == 2 || Acase == 3)
			{
				if (IS_VALID_backward(allMVCPoints, CLC, phi, _nearInterval, allMVCPoints.front()(1)))
					endInterval(1) = allMVCPoints.front()(1);
				else
				{
					LOG("Case B3b, isValid is false, the path is not traversable. AVP backward failure.");
					return false;
				}
			}
		}
		else if (Bcase == 4)
		{
			if (IS_VALID_backward(allMVCPoints, CLC, phi, _nearInterval, allMVCPoints.front()(1)))
				endInterval(1) = allMVCPoints.front()(1);
			else
			{
				LOG("Case B4, isValid is false, the path is not traversable. AVP backward failure.");
				return false;
			}
		}

		/* Step C : determining the minimum final velocity */
		endInterval(0) = findsdotminBybinearSearch(allMVCPoints, CLC, phi, _nearInterval, endInterval(1), BACKWARD);
		std::cout << endInterval << std::endl;

		savedata(allMVCPoints, allSwitchPoint, LC_copy, CLC, phi);

		return true;
	}

	void AVP_RRT::settingtopp(std::list<VectorX>& Pnew)
	{
		_topp->initialization();
		unsigned int RowSize = Pnew.front().size();
		unsigned int ColSize = Pnew.size();
		MatrixX q_data(RowSize, ColSize);
		unsigned int cnt = 0;
		for (std::list<VectorX>::iterator it = Pnew.begin(); it != Pnew.end(); it++)
			q_data.col(cnt++) = (*it);

		_topp->setJointTrajectory(q_data);
	}

	bool AVP_RRT::calculateLimitingCurves(const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints,
		const std::vector<unsigned int>& allMVCPointsFlag, const std::vector<SwitchPoint>& allSwitchPoint,
		std::vector<std::list<Vector2, Eigen::aligned_allocator<Vector2>>>& LC)
	{
		Real s_cur, sdot_cur, s_swi, sdot_swi;
		Real ds = _topp->getStepSize();
		Real si = _topp->getInitialParam();
		Real sf = _topp->getFinalParam();
		unsigned int numOfswi = allSwitchPoint.size();
		int swi;

		std::list<Vector2, Eigen::aligned_allocator<Vector2>> tmp_list;

		for (unsigned int i = 0; i < numOfswi; i++)
		{
			std::cout << i << std::endl;


			tmp_list.clear();
			LC.push_back(tmp_list);
			s_swi = round(allSwitchPoint[i]._s / ds) * ds;
			sdot_swi = allSwitchPoint[i]._sdot;

			// backward integration
			s_cur = s_swi; sdot_cur = sdot_swi;
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
			backward_fcnpointer = &AVP_RRT::backwardInt;
			while (s_cur > si)
			{
				swi = (this->*backward_fcnpointer)(s_cur, sdot_cur, LC[i], allMVCPoints, allMVCPointsFlag);
				if (swi == 0)
					return false;
				else if (swi == 1)
					break;
			}

			// forward integration
			s_cur = s_swi; sdot_cur = sdot_swi;
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
			forward_fcnpointer = &AVP_RRT::forwardInt;
			while (s_cur < sf)
			{
				swi = (this->*forward_fcnpointer)(s_cur, sdot_cur, LC[i], allMVCPoints, allMVCPointsFlag);
				if (swi == 0)
					return false;
				else if (swi == 1)
					break;
			}
		}
		return true;
	}

	int AVP_RRT::forwardInt(Real& s_cur, Real& sdot_cur, std::list<Vector2, Eigen::aligned_allocator<Vector2>>& LC,
		const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints, const std::vector<unsigned int>& allMVCPointsFlag)
	{
		Real beta_cur, sdot_MVC;
		unsigned int idx;
		int flag;

		beta_cur = _topp->determineAlphaBeta(s_cur, sdot_cur)(1);
		_topp->forwardIntegrate(s_cur, sdot_cur, beta_cur);
		if (s_cur > sf)
			return 1;
		idx = round(s_cur / ds);
		sdot_MVC = allMVCPoints[idx](1);
		flag = allMVCPointsFlag[idx];

		if (sdot_cur < 1e-5)
			return 0;
		else if (sdot_cur > sdot_MVC)
		{
			if (flag == 2)
				return 1;
			else if (flag == 1)
				forward_fcnpointer = &AVP_RRT::forwardIntVel;
		}
		
		LC.push_back(Vector2(s_cur, sdot_cur));

		return (this->*forward_fcnpointer)(s_cur, sdot_cur, LC, allMVCPoints, allMVCPointsFlag);
	}

	int AVP_RRT::forwardIntVel(Real& s_cur, Real& sdot_cur, std::list<Vector2, Eigen::aligned_allocator<Vector2>>& LC,
		const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints, const std::vector<unsigned int>& allMVCPointsFlag)
	{
		Vector2 alphabeta;
		Real sdot_MVC, s_next, sdot_next, ds, slope;
		Real alpha_cur, beta_cur;
		unsigned int idx;
		int flag;

		//alphabeta = _topp->determineAlphaBeta(s_cur, sdot_cur);
		//alpha_cur = alphabeta(0);
		//beta_cur = alphabeta(1);

		s_next = s_cur + ds;
		if(s_next > sf)
			return 1;

		idx = round(s_next / ds);
		sdot_next = allMVCPoints[idx](1);
		flag = allMVCPointsFlag[idx];
		if (flag == 2)
			return 1;

		//slope = (sdot_next - sdot_cur) / (s_next - s_cur);


		//if (beta_cur >= slope && slope >= alpha_cur)
		//{
		//	s_cur = s_next; sdot_cur = sdot_next;
		//	LC.push_back(Vector2(s_cur, sdot_cur));
		//}
		//else if (beta_cur < slope)
		//	forward_fcnpointer = &AVP_RRT::forwardInt;
		//if (beta_cur >= slope && slope >= alpha_cur)
		//{
		//	s_cur = s_next; sdot_cur = sdot_next;
		//	LC[i].push_back(Vector2(s_cur, sdot_cur));
		//}
		//else if (beta_cur < slope)
		//{
		//	// forward integrate..,.
		//}

		LC.push_back(Vector2(s_cur, sdot_cur));
		return (this->*forward_fcnpointer)(s_cur, sdot_cur, LC, allMVCPoints, allMVCPointsFlag);
	}

	int AVP_RRT::backwardInt(Real& s_cur, Real& sdot_cur, std::list<Vector2, Eigen::aligned_allocator<Vector2>>& LC,
		const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints, const std::vector<unsigned int>& allMVCPointsFlag)
	{
		Real alpha_cur, sdot_MVC;
		unsigned int idx;
		int flag;

		alpha_cur = _topp->determineAlphaBeta(s_cur, sdot_cur)(0);
		_topp->backwardIntegrate(s_cur, sdot_cur, alpha_cur);
		if (s_cur < _topp->getInitialParam())
		//if (s_cur < si)
			return 1;
		idx = round(s_cur / _topp->getStepSize());
		sdot_MVC = allMVCPoints[idx](1);
		flag = allMVCPointsFlag[idx];

		if (sdot_cur < 1e-5)
			return 0;
		else if (sdot_cur > sdot_MVC)
		{
			if (flag == 2)
				return 1;
			else if (flag == 1)
				backward_fcnpointer = &AVP_RRT::backwardIntVel;
		}

		LC.push_front(Vector2(s_cur, sdot_cur));
		return (this->*backward_fcnpointer)(s_cur, sdot_cur, LC, allMVCPoints, allMVCPointsFlag);
	}
	int AVP_RRT::backwardIntVel(Real& s_cur, Real& sdot_cur, std::list<Vector2, Eigen::aligned_allocator<Vector2>>& LC,
		const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints, const std::vector<unsigned int>& allMVCPointsFlag)
	{
		Vector2 alphabeta;
		Real sdot_MVC, s_next, sdot_next, ds, slope;
		Real alpha_cur, beta_cur;
		unsigned int idx;
		int flag;

		ds = _topp->getStepSize();

		//alphabeta = _topp->determineAlphaBeta(s_cur, sdot_cur);
		//alpha_cur = alphabeta(0);
		//beta_cur = alphabeta(1);

		s_next = s_cur - ds;
		if (s_next < _topp->getInitialParam())
			return 1;

		idx = round(s_next / ds);
		sdot_next = allMVCPoints[idx](1);
		flag = allMVCPointsFlag[idx];
		if (flag == 2)
			return 1;
		
		//slope = (sdot_next - sdot_cur) / (s_next - s_cur);

		s_cur = s_next; sdot_cur = sdot_next;

		//if (beta_cur >= slope && slope >= alpha_cur)
		//{
		//	s_cur = s_next; sdot_cur = sdot_next;
		//	LC[i].push_back(Vector2(s_cur, sdot_cur));
		//}
		//else if (beta_cur < slope)
		//{
		//	// forward integrate..,.
		//}

		LC.push_front(Vector2(s_cur, sdot_cur));
		return (this->*backward_fcnpointer)(s_cur, sdot_cur, LC, allMVCPoints, allMVCPointsFlag);
	}

	void AVP_RRT::calculateCLC(std::vector<std::list<Vector2, Eigen::aligned_allocator<Vector2>>>& LC,
		std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& CLC)
	{
		// 0.199 ���� �̻���..

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

		for (unsigned int i = 0; i < LCsize; i++)
		{
			std::list<Vector2, Eigen::aligned_allocator<Vector2>>::iterator it_tmp;
			it_tmp = LC[i].begin();
			it_begin.push_back(it_tmp);
		}

		// find s_end
		Real tmp = -std::numeric_limits<Real>::max();
		for (unsigned int i = 0; i < LCsize; i++)
		{
			if ((*(it_end[i]))(0) > tmp)
				tmp = (*(it_end[i]))(0);
		}
		s_end = tmp;

		// find s_begin
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
		Real ds = _topp->getStepSize();

		Real tmp_min;
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
					idx.push_back(i);
			}

			tmp_min = std::numeric_limits<Real>::max();
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

		idx.clear();

		s_cur = s_end;

		// for last value
		for (unsigned int i = 0; i < LCsize; i++)
		{
			if ((*(it_end[i]))(0) < s_end)
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

		tmp_min = std::numeric_limits<Real>::max();
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

	unsigned int AVP_RRT::determineAresult(const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints,
		const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& CLC, Real& sdot_beg_star, AVPFLAG avpflag)
	{
		bool si_col = false, sf_col = false;

		Real eps_init = std::pow((CLC.front()(0) - _topp->getInitialParam()),2);
		Real eps_final = std::pow((CLC.back()(0) - _topp->getFinalParam()), 2);

		if (eps_init <= 1e-5)
			si_col = true;
		if (eps_final <= 1e-5)
			sf_col = true;

		if (si_col == false && sf_col == false)
		{
			if (avpflag == FORWARD)
				sdot_beg_star = allMVCPoints.front()(1);
			else if (avpflag == BACKWARD)
				sdot_beg_star = allMVCPoints.back()(1);
			return 2;
		}
		else if (si_col == true && sf_col == false)
		{
			if (avpflag == FORWARD)
			{
				sdot_beg_star = CLC.front()(1);
				return 3;
			}
			else if (avpflag == BACKWARD)
			{
				sdot_beg_star = allMVCPoints.back()(1);
				return 4;
			}
		}
		else if (si_col == false && sf_col == true)
		{
			if (avpflag == FORWARD)
			{
				sdot_beg_star = allMVCPoints.front()(1);
				return 4;
			}
			else if (avpflag == BACKWARD)
			{
				sdot_beg_star = CLC.back()(1);
				return 3;
			}
		}
		else if (si_col == true && sf_col == true)
		{
			if (avpflag == FORWARD)
				sdot_beg_star = CLC.front()(1);
			else if (avpflag == BACKWARD)
				sdot_beg_star = CLC.back()(1);
			return 5;
		}
	}

	unsigned int AVP_RRT::determineAVPBresult(const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints,
		const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& CLC, const Real sdot_init, 
		std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& phi)
	{
		Real s_cur, sdot_cur, beta_cur, sf, ds, sdot_MVC, sdot_CLC;
		int idx, flag;
		sf = _topp->getFinalParam();
		ds = _topp->getStepSize();
		s_cur = _topp->getInitialParam();
		sdot_cur = sdot_init;
		phi.push_back(Vector2(s_cur, sdot_cur));

		while (true)
		{
			beta_cur = _topp->determineAlphaBeta(s_cur, sdot_cur)(1);
			_topp->forwardIntegrate(s_cur, sdot_cur, beta_cur);
			idx = round(s_cur / ds);
			sdot_MVC = allMVCPoints[idx](1);
			sdot_CLC = CLC[idx](1);

			if (sdot_cur < 0)
				return 1;
			else if (s_cur > sf)
				return 2;
			else if (sdot_cur > sdot_CLC)
				return 3;
			else if (sdot_cur > sdot_MVC)
				return 4;

			phi.push_back(Vector2(s_cur, sdot_cur));
		}
	}

	unsigned int AVP_RRT::determineAVPBackwardBresult(const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints,
		const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& CLC, const Real sdot_init,
		std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& phi)
	{
		std::list<Vector2, Eigen::aligned_allocator<Vector2>> phi_tmp;
		Real s_cur, sdot_cur, alpha_cur, si, ds, sdot_MVC, sdot_CLC;
		int idx, flag;
		si = _topp->getInitialParam();
		ds = _topp->getStepSize();
		s_cur = _topp->getFinalParam() - 1e-6;
		sdot_cur = sdot_init;
		phi_tmp.push_front(Vector2(s_cur, sdot_cur));

		unsigned int Bresult;

		while (true)
		{
			alpha_cur = _topp->determineAlphaBeta(s_cur, sdot_cur)(0);
			_topp->backwardIntegrate(s_cur, sdot_cur, alpha_cur);
			idx = round(s_cur / ds);
			sdot_MVC = allMVCPoints[idx](1);
			sdot_CLC = CLC[idx](1);

			if (sdot_cur < 0)
			{
				Bresult = 1;
				break;
			}
			else if (s_cur < si)
			{
				Bresult = 2;
				break;
			}

			else if (sdot_cur > sdot_CLC)
			{
				Bresult = 3;
				break;
			}
			else if (sdot_cur > sdot_MVC)
			{
				Bresult = 4;
				break;
			}
			phi_tmp.push_front(Vector2(s_cur, sdot_cur));
		}

		for (std::list<Vector2, Eigen::aligned_allocator<Vector2>>::iterator it = phi_tmp.begin(); it != phi_tmp.end(); it++)
			phi.push_back(*it);

		return Bresult;
	}

	bool AVP_RRT::IS_VALID(const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints,
		const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& CLC, const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& phi, 
		const Vector2& nearInterval, const Real sdot_test)
	{
		Real s_cur, sdot_cur, alpha_cur, si, ds, sdot_phi, sdot_CLC;
		int idx, flag;
		si = _topp->getInitialParam();
		ds = _topp->getStepSize();
		s_cur = _topp->getFinalParam();
		sdot_cur = sdot_test;

		int phiSize = phi.size();

		while (true)
		{
			alpha_cur = _topp->determineAlphaBeta(s_cur, sdot_cur)(0);
			_topp->backwardIntegrate(s_cur, sdot_cur, alpha_cur);
			idx = round(s_cur / ds);
			sdot_CLC = CLC[idx](1);
			if (idx > phiSize)
				sdot_phi = std::numeric_limits<Real>::max();
			else
				sdot_phi = phi[idx](1);

			if (sdot_cur < 0)
				return false;
			else if (s_cur < si)
			{
				if (sdot_cur < nearInterval(0))
					return false;
				else if (sdot_cur > nearInterval(0))
					return true;
			}
			else if (sdot_cur > sdot_phi)
				return true;
			else if (sdot_cur > sdot_CLC)
				return true;
		}
		return false;
	}

	bool AVP_RRT::IS_VALID_backward(const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints,
		const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& CLC, const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& phi,
		const Vector2& nearInterval, const Real sdot_test)
	{
		Real s_cur, sdot_cur, beta_cur, sf, ds, sdot_phi, sdot_CLC;
		int idx, flag;
		sf = _topp->getFinalParam();
		ds = _topp->getStepSize();
		s_cur = _topp->getInitialParam();
		sdot_cur = sdot_test;

		Real s_phi_init = phi.front()(0);
		unsigned int cnt = 0;

		while (true)
		{
			beta_cur = _topp->determineAlphaBeta(s_cur, sdot_cur)(1);
			_topp->forwardIntegrate(s_cur, sdot_cur, beta_cur);
			idx = round(s_cur / ds);
			sdot_CLC = CLC[idx](1);
			if (s_cur < s_phi_init)
				sdot_phi = std::numeric_limits<Real>::max();
			else if (s_cur >= s_phi_init)
				sdot_phi = phi[cnt++](1);
				

			if (sdot_cur < 0)
				return false;
			else if (s_cur > sf)
			{
				if (sdot_cur < nearInterval(0))
					return false;
				else if (sdot_cur > nearInterval(0))
					return true;
			}
			else if (sdot_cur > sdot_phi)
				return true;
			else if (sdot_cur > sdot_CLC)
				return true;
		}
		return false;
	}

	Real AVP_RRT::findsdotminBybinearSearch(const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints,
		const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& CLC, const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& phi,
		const Vector2& nearInterval, const Real sdot_test, AVPFLAG avpflag)
	{
		Real sdot_max = sdot_test;
		Real sdot_min = 0;

		Real sdot_before;
		Real sdot_up = sdot_max;
		Real sdot_down = sdot_min;
		Real sdot_cur = sdot_min;
		Real epsilon = sdot_max * 1e-3;
		unsigned int max_iter = 100;
		unsigned int cnt = 0;
		bool isValid;

		if (avpflag == FORWARD)
			isValid = IS_VALID(allMVCPoints, CLC, phi, nearInterval, sdot_cur);
		else if (avpflag == BACKWARD)
			isValid = IS_VALID_backward(allMVCPoints, CLC, phi, nearInterval, sdot_cur);
		if (isValid)
			return 0;

		while (true)
		{
			if (isValid)
			{
				sdot_before = sdot_cur;
				sdot_up = sdot_cur;
				sdot_cur = (sdot_cur + sdot_down) * 0.5;
			}
			else
			{
				sdot_before = sdot_cur;
				sdot_down = sdot_cur;
				sdot_cur = (sdot_cur + sdot_up) * 0.5;
			}

			if (avpflag == FORWARD)
				isValid = IS_VALID(allMVCPoints, CLC, phi, nearInterval, sdot_cur);
			else if (avpflag == BACKWARD)
				isValid = IS_VALID_backward(allMVCPoints, CLC, phi, nearInterval, sdot_cur);

			if (isValid == true && cnt > max_iter)
				break;
			else if (isValid == true && std::abs(sdot_before - sdot_cur) < epsilon)
				break;

			cnt++;

		}
		return sdot_cur;
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

		return _nodes[index]; // nearest node�� �´°� -> ������ϸ鼭 Ȯ���غ���.
	}

	// addVertex �θ��� ���� vertex ������� �� �� �����Ѵ����� �ֱ�

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