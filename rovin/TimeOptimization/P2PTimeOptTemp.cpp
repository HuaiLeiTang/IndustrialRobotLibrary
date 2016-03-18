#include "P2PTimeOptTemp.h"

using namespace std;

namespace rovin
{
	AVP_RRT::AVP_RRT(const SerialOpenChainPtr& robot, CONSTRAINT_TYPE constraintType) : _robot(robot), _constraintType(constraintType)
	{
		_dof = robot->getNumOfJoint();
		Real ds = 0.5e-2, vi = 0, vf = 0, si = 0, sf = 1;
		_topp = TOPPPtr(new TOPP(_robot, vi, vf, ds, si, sf, constraintType));

		_ds = _topp->getStepSize();
		_si = _topp->getInitialParam();
		_sf = _topp->getFinalParam();

		_dof = robot->getNumOfJoint();
		//_stepsize = 0.1;
		_stepsize = 0.5;
		_curSegment = 0;

		srand((int)time(NULL));
	}

	void AVP_RRT::generateTrajectory()
	{
		// initialization! -> way point setting
		_numSegment = _waypoints.size() - 1;
		_segmentPath.resize(_numSegment);

		_curSegment = 0;
		RETURNFLAG retFlag;
		do {
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
			std::cout << iter << std::endl;
			makeRandomConfig(qrand);
			//qrand(0) = 2.71335; qrand(1) = 1.10668;  qrand(2) = -1.64206;
			//qrand(3) = -0.280439; qrand(4) = -1.94905;  qrand(5) = -4.18502;


			std::cout << "qrand : " << qrand << std::endl;

			//if collisionchekc(qrand)
			//	continue;


			// extend from start tree to random configuration
			candiVertex = new Vertex();
			extendTree(&_startTree, qrand, true, &candiVertex);
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
			candiVertex = new Vertex();
			extendTree(&_goalTree, qrand, false, &candiVertex);
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

		} while (++iter< AVP_RRT_MAX_ITER);

		if (iter >= AVP_RRT_MAX_ITER)
			flag = EXCEED_MAX_ITER;


		//if (flag == SUCCESS)
		//{
		//	_wayPointInterval = _goalTree._nodes[0]->_interval;
		//}


		_startTree.clearTree();
		_goalTree.clearTree();
		return flag;

	}


	void AVP_RRT::treeInitialization(int idx)
	{
		Vector2 tmp(0, std::numeric_limits<Real>::max());

		Vertex * startVertex = new Vertex(_waypoints[idx].getJointq(), _waypoints[idx].getJointqdot(), tmp, std::list<VectorX>(), NULL);
		Vertex * goalVertex = new Vertex(_waypoints[idx + 1].getJointq(), _waypoints[idx + 1].getJointqdot(), tmp, std::list<VectorX>(), NULL);
		
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
			qrand(i) = (Real)(((ub - lb)*rand()) / (RAND_MAX + 0.0) + lb);
		}
	}

	bool AVP_RRT::extendTree(Tree * tree, const VectorX qrand, bool forward, Vertex ** candiVertex) // 반환 Vertex* 로??
	{
		// find nearest vertex
		double tmpDist;
		Vertex * nVertex = tree->findNearestNeighbor(qrand, &tmpDist);

		std::cout << "nVertex config : " << nVertex->_config << std::endl;




		// interpolate between nVertex.config and qnew(configuration far away from nVertex.config about stepsize)
		std::list<VectorX> Pnew;
		VectorX qnew, qvel;
		if (!interpolate(nVertex, qrand, tmpDist, forward, Pnew, qnew, qvel))
		{
			//FAILURE!!
		}




		// run AVP algorithm
		Vector2 endInterval;
		bool isSucceeded;
		if (forward)
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
		{
			(*candiVertex) = new Vertex(qnew, qvel, endInterval, Pnew, nVertex);
			return true;
		}
		else
		{
			(*candiVertex) = NULL;
			return false;
		}

	}

	bool AVP_RRT::interpolate_tmp(Vertex * nVertex, const VectorX qrand, const Real dist, bool forward, std::list<VectorX> & Pnew, VectorX & qnew, VectorX & qvel)
	{
		VectorX qnear = nVertex->_config;
		VectorX qnearvel = nVertex->_configVel;

		qnew = qrand;
		//if (dist < _stepsize)
		//	qnew = qrand;
		//else
		//	qnew = qnear + (qrand - qnear)*_stepsize / dist;

		cout << "qnew : " << qnew << endl;

		unsigned int numOfCP = 6;
		unsigned int degree = 3;
		unsigned int order = degree + 1;
		unsigned int numOfknot = numOfCP + degree + 1;

		VectorX knot(numOfknot);
		MatrixX CP(_dof, numOfCP);

		// make knot
		for (unsigned int i = 0; i < order; i++)
		{
			knot[i] = _si;
			knot[numOfknot - i - 1] = _sf;
		}
		Real delta = (_sf - _si) / (numOfCP - degree);
		for (unsigned int i = 0; i < (numOfCP - degree); i++)
			knot[order + i] = delta * (i + 1);

		// make boundary condition
		CP.col(0) = qnear;
		CP.col(1) = delta / degree * qnearvel + CP.col(0);
		CP.col(CP.cols() - 1) = qnew;

		// make control point
		VectorX delta_cp = (CP.col(CP.cols() - 1) - CP.col(1)) / double((numOfCP - 2));
		for (int i = 2; i < CP.cols() - 1; i++)
			CP.col(i) = CP.col(i - 1) + delta_cp;

		BSpline<-1, -1, -1> q(knot, CP);
		BSpline<-1, -1, -1> qs = q.derivative();

		unsigned int numOfdata = 500;
		Real delta_s = (_sf - _si) / double(numOfdata - 1);
		Real s_cur = _si;

		std::vector<Real> ss;

		Pnew.push_back(q(_si));
		for (unsigned int i = 1; i < numOfdata - 1; i++)
		{
			s_cur += delta_s;
			Pnew.push_back(q(s_cur));
		}
		Pnew.push_back(q(_sf - 1e-5));


		qvel = qs(_sf - 1e-5);

		std::cout << "qvel : " << qvel << std::endl;


		return false;

	}


	bool AVP_RRT::interpolate(Vertex * nVertex, const VectorX qrand, const Real dist, bool forward, std::list<VectorX> & Pnew, VectorX & qnew, VectorX & qvel)
	{
		VectorX qnear = nVertex->_config;
		VectorX qnearvel = nVertex->_configVel;

		//qnew = qrand;
		if (dist < _stepsize)
			qnew = qrand;
		else
			qnew = qnear + (qrand - qnear)*_stepsize / dist;

		cout << "qnew : " << qnew << endl;

		unsigned int numOfCP = 6;
		unsigned int degree = 3;
		unsigned int order = degree + 1;
		unsigned int numOfknot = numOfCP + degree + 1;

		VectorX knot(numOfknot);
		MatrixX CP(_dof, numOfCP);

		// make knot
		for (unsigned int i = 0; i < order; i++)
		{
			knot[i] = _si;
			knot[numOfknot - i - 1] = _sf;
		}
		Real delta = (_sf - _si) / (numOfCP - degree);
		for (unsigned int i = 0; i < (numOfCP - degree); i++)
			knot[order + i] = delta * (i + 1);

		// make boundary condition
		CP.col(0) = qnear;
		CP.col(1) = delta / degree * qnearvel + CP.col(0);
		CP.col(CP.cols() - 1) = qnew;

		// make control point
		VectorX delta_cp = (CP.col(CP.cols() - 1) - CP.col(1)) / double((numOfCP - 2));
		for (int i = 2; i < CP.cols() - 1; i++)
			CP.col(i) = CP.col(i - 1) + delta_cp;

		BSpline<-1, -1, -1> q(knot, CP);
		BSpline<-1, -1, -1> qs = q.derivative();

		unsigned int numOfdata = 500;
		Real delta_s = (_sf - _si) / double(numOfdata - 1);
		Real s_cur = _si;

		std::vector<Real> ss;

		Pnew.push_back(q(_si));
		for (unsigned int i = 1; i < numOfdata - 1; i++)
		{
			s_cur += delta_s;
			Pnew.push_back(q(s_cur));
		}
		Pnew.push_back(q(_sf - 1e-5));

		//for (std::list<VectorX>::iterator it = Pnew.begin(); it != Pnew.end(); it++)
		//	std::cout << "Pnew" << *it << std::endl;

		qvel = qs(_sf - 1e-5);

		std::cout << "qvel : " << qvel << std::endl;


		// feasibility test

		bool testCollision = true; // true: succeded (no collision), false: collision detected
								   // collision test routine for Pnew is needed...


								   // if nearest vertex is the root vertex
		//VectorX qs; // calc qs is needed
		VectorX qs_zero = qs(_si);

		std::cout << "qs_zero : " << qs_zero << std::endl;

		bool testRoot;
		if (nVertex->_parentVertex == NULL)
		{
			if (forward)
				testRoot = testRootVertex(Pnew, qs_zero);
			else
				testRoot = testRootVertexbackward(Pnew, qs_zero);
		}

		if (testCollision && testRoot)
			return true;

		return false;

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
			VectorX qtmp, qvel;
			if (!interpolate(vertex, (*oVertex)->_config, tmpDist, forward, Ptmp, qtmp, qvel))
				return false;

			Vector2 endInterval;
			if(forward) // vertex belongs to start tree, and try to connect to goal tree
			{
				if (runAVP(Ptmp, vertex->_interval, endInterval))
				{
					if (checkIntersectionOfTwoIntervals(endInterval, (*oVertex)->_interval))
					{
						(*cVertex)->_config = (*oVertex)->_config;
						(*cVertex)->_configVel = qvel;
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
						(*cVertex)->_configVel = qvel;
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



	bool AVP_RRT::testRootVertex(std::list<VectorX>& Pnew, VectorX& qs)
	{
		settingtopp(Pnew);

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

		saveMVC(allMVCPoints);
		saveSwitchPoint(allSwitchPoint);
		saveIntVector2txt(allMVCPointsFlag, "C:/Users/crazy/Desktop/Time optimization/avp test/MVC_flag.txt");
		
		std::cout << allMVCPoints.size() << std::endl;
		std::cout << allMVCPointsFlag.size() << std::endl;
		std::cout << allSwitchPoint.size() << std::endl;

		A1switch = calculateLimitingCurves(allMVCPoints, allMVCPointsFlag, allSwitchPoint, LC);

		std::vector<std::list<Vector2, Eigen::aligned_allocator<Vector2>>> LC_copy;
		LC_copy = LC; ///< for save LC

		if (A1switch == false)
		{
			LOG("A1 case, Limiting curve reaches sdot = 0. AVP failure.");
			return false;
		}

		calculateCLC(LC, CLC);

		std::cout << "CLC size : " << CLC.size() << std::endl;

		unsigned int Acase;
		Real sdot_beg_star;
		Acase = determineAresult(allMVCPoints, CLC, sdot_beg_star, FORWARD);

		
		// feasibility test of sdot_beg_star
		VectorX tmpVec(_dof);
		tmpVec = _waypoints[_curSegment].getJointqdot();
		tmpVec = tmpVec.cwiseProduct(qs.cwiseInverse());

		if (sdot_beg_star < tmpVec.maxCoeff())
			return false; // failure case
		
		return true;

	}
	bool AVP_RRT::testRootVertexbackward(std::list<VectorX>& Pnew, VectorX& qs)
	{
		settingtopp(Pnew);

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

		if (A1switch == false)
		{
			LOG("A1 case, Limiting curve reaches sdot = 0. AVP failure.");
			return false;
		}

		calculateCLC(LC, CLC);
		unsigned int Acase;
		Real sdot_end_star;
		Acase = determineAresult(allMVCPoints, CLC, sdot_end_star, BACKWARD);

		// feasibility test of sdot_beg_star
		VectorX tmpVec(_dof);
		tmpVec = _waypoints[_curSegment+1].getJointqdot();
		tmpVec.cwiseProduct(qs.cwiseInverse());

		if (sdot_end_star < tmpVec.maxCoeff())
			return false; // failure case

		return true;
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
		bool A1switch = true;

		_topp->calculateAllMVCPoint();
		_topp->calculateAllSwitchPoint();
		allMVCPoints = _topp->getAllMVCPoint();
		allMVCPointsFlag = _topp->getAllMVCPointFlag();
		allSwitchPoint = _topp->getAllSwitchPoint();
		
		if (allSwitchPoint.size() != 0)
			A1switch = calculateLimitingCurves(allMVCPoints, allMVCPointsFlag, allSwitchPoint, LC);

		std::vector<std::list<Vector2, Eigen::aligned_allocator<Vector2>>> LC_copy;
		LC_copy = LC; ///< for save LC

		if (A1switch == false)
		{
			LOG("A1 case, Limiting curve reaches sdot = 0. AVP failure.");
			return false;
		}

		if (allSwitchPoint.size() != 0)
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

		saveData(allMVCPoints, allSwitchPoint, LC_copy, CLC, phi);

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

		Real sdot_end_max_star = std::min(sdot_end_max, sdot_end_star);  ///< 이게 맞는건지 모르겠음..
		Vector2 _nearInterval(sdot_end_min, sdot_end_max_star); ///< 이게 맞는건지 모르겠음..
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

		saveData(allMVCPoints, allSwitchPoint, LC_copy, CLC, phi);
		
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
		unsigned int numOfswi = allSwitchPoint.size();
		unsigned int swi;

		Vector2 alphabeta;
		std::list<Vector2, Eigen::aligned_allocator<Vector2>> backward_list;
		std::list<Vector2, Eigen::aligned_allocator<Vector2>> forward_list;

		for (unsigned int i = 0; i < numOfswi; i++)
		{
			backward_list.clear();
			forward_list.clear();
			s_swi = round(allSwitchPoint[i]._s / _ds) * _ds;
			sdot_swi = allSwitchPoint[i]._sdot;

			// backward integration
			s_cur = s_swi;
			sdot_cur = sdot_swi;
			backward_list.push_front(Vector2(s_cur, sdot_cur));

			std::cout << "switch point id : " << allSwitchPoint[i]._id << std::endl;

			if (allSwitchPoint[i]._id == SwitchPoint::SINGULAR)
			{
				Real lambda = allSwitchPoint[i]._lambda;
				unsigned int numOfSPInt = 3;
				for (unsigned int j = 0; j < numOfSPInt; j++)
				{
					s_cur -= _ds;
					sdot_cur -= lambda * _ds;
					backward_list.push_front(Vector2(s_cur, sdot_cur));
				}
			}

			swi = backwardIntegrate(s_cur, sdot_cur, backward_list, allMVCPoints, allMVCPointsFlag);
			
			if (swi == 0)
			{
				s_swi = backward_list.back()(0);
				sdot_swi = backward_list.back()(1);
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
					s_cur += _ds;
					sdot_cur += lambda * _ds;
					forward_list.push_back(Vector2(s_cur, sdot_cur));
				}
			}

			while (true)
			{
				swi = forwardIntegrate(s_cur, sdot_cur, forward_list, allMVCPoints, allMVCPointsFlag);
			
				if (swi == 1)
					break;
				else if (swi == 0)
				{
					backward_list.clear();
					s_cur = forward_list.front()(0);
					sdot_cur = forward_list.front()(1);
				}
				
				swi = backwardIntegrate(s_cur, sdot_cur, backward_list, allMVCPoints, allMVCPointsFlag);
				
				if (swi == 1)
					break;
				else if (swi == 0)
				{
					forward_list.clear();
					s_cur = backward_list.back()(0);
					sdot_cur = backward_list.back()(1);
				}
			}

			backward_list.insert(backward_list.end(), forward_list.begin(), forward_list.end());
			LC.push_back(backward_list);
		}
		return true;
	}

	unsigned int AVP_RRT::forwardIntegrate(Real& s_cur, Real& sdot_cur, std::list<Vector2, Eigen::aligned_allocator<Vector2>>& LC,
		const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints,
		const std::vector<unsigned int>& allMVCPointsFlag)
	{
		unsigned int swi_forward;
		unsigned int swi_forwardVel = 2;
		unsigned int swi_forwardback = 1;
		while (true)
		{
			if (swi_forwardVel == 2)
				swi_forward = forwardInt(s_cur, sdot_cur, LC, allMVCPoints, allMVCPointsFlag);

			if (swi_forward == 0)
				return false;
			else if (swi_forward == 1)
				break;
			else if (swi_forward == 2)
				swi_forwardVel = forwardIntVel(s_cur, sdot_cur, LC, allMVCPoints, allMVCPointsFlag);

			if (swi_forwardVel == 1)
				break;
			else if (swi_forwardVel == 3)
			{
				swi_forwardback = forwardbackInt(s_cur, sdot_cur, LC, allMVCPoints, allMVCPointsFlag);
				swi_forward = 2;
			}
		}
		return swi_forwardback;
	}

	unsigned int AVP_RRT::forwardInt(Real& s_cur, Real& sdot_cur, std::list<Vector2, Eigen::aligned_allocator<Vector2>>& LC, 
		const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints,
		const std::vector<unsigned int>& allMVCPointsFlag)
	{
		Real beta_cur, sdot_MVC;
		unsigned int idx;
		int flag;
		while (s_cur < _sf)
		{
			beta_cur = _topp->determineAlphaBeta(s_cur, sdot_cur)(1);
			_topp->forwardIntegrate(s_cur, sdot_cur, beta_cur);
			if (sdot_cur < 1e-4)
				return 0;

			idx = (unsigned int)round(s_cur / _ds);
			sdot_MVC = allMVCPoints[idx](1);
			flag = allMVCPointsFlag[idx];

			if (sdot_cur > sdot_MVC)
			{
				sdot_cur = sdot_MVC;
				LC.push_back(Vector2(s_cur, sdot_cur));
				if (flag == 2)
					return 1;
				if (flag == 1)
					return 2;
			}
			LC.push_back(Vector2(s_cur, sdot_cur));
		}
		return 1;
	}

	unsigned int AVP_RRT::forwardIntVel(Real& s_cur, Real& sdot_cur, std::list<Vector2, Eigen::aligned_allocator<Vector2>>& LC,
		const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints,
		const std::vector<unsigned int>& allMVCPointsFlag)
	{
		Vector2 alphabeta;
		Real s_next, sdot_next, alpha_cur, beta_cur, slope;
		unsigned int idx;
		int flag;

		while (true)
		{
			s_next = s_cur + _ds;
			if (s_next > _sf)
			{
				s_cur = _sf; sdot_cur = allMVCPoints[(unsigned int)round(_sf / _ds)](1);
				LC.push_back(Vector2(s_cur, sdot_cur));
				return 1;
			}

			idx = (unsigned int)round(s_next / _ds);
			sdot_next = allMVCPoints[idx](1);
			flag = allMVCPointsFlag[idx];
			if (flag == 2)
				return 1;

			alphabeta = _topp->determineAlphaBeta(s_cur, sdot_cur);
			alpha_cur = alphabeta(0);
			beta_cur = alphabeta(1);

			slope = (sdot_next - sdot_cur) / (s_next - s_cur);

			if (beta_cur >= slope && slope >= alpha_cur)
			{
				s_cur = s_next; sdot_cur = sdot_next;
				LC.push_back(Vector2(s_cur, sdot_cur));
			}
			else if (beta_cur < slope) 
				return 2;
			else if (slope < alpha_cur)
				return 3;
		}
	}

	unsigned int AVP_RRT::forwardbackInt(Real& s_cur, Real& sdot_cur, std::list<Vector2, Eigen::aligned_allocator<Vector2>>& LC,
		const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints,
		const std::vector<unsigned int>& allMVCPointsFlag)
	{
		Vector2 alphabeta;
		Real s_next, sdot_next, alpha_cur, beta_cur, slope;
		unsigned int idx;
		int flag;

		while (true)
		{
			s_next = s_cur + _ds;
			if (s_next > _sf)
			{
				s_cur = _sf; sdot_cur = allMVCPoints[(unsigned int)round(_sf / _ds)](1);
				break;
			}
			idx = (unsigned int)round(s_next / _ds);
			sdot_next = allMVCPoints[idx](1);
			flag = allMVCPointsFlag[idx];
			if (flag == 2)
				break;
			alphabeta = _topp->determineAlphaBeta(s_cur, sdot_cur);
			alpha_cur = alphabeta(0);
			beta_cur = alphabeta(1);
			slope = (sdot_next - sdot_cur) / (s_next - s_cur);
			if (beta_cur >= slope && slope >= alpha_cur)
				break;
			s_cur = s_next; sdot_cur = sdot_next;
		}
		
		// backward integration
		unsigned int result;
		Real sdot_LC;
		Real s_LC_back = LC.back()(0);
		std::list<Vector2, Eigen::aligned_allocator<Vector2>> tmp_back_list;
		tmp_back_list.push_front(Vector2(s_cur, sdot_cur));
		while (true)
		{
			alpha_cur = _topp->determineAlphaBeta(s_cur, sdot_cur)(0);
			_topp->backwardIntegrate(s_cur, sdot_cur, alpha_cur);
			
			if (s_LC_back >= s_cur)
			{
				sdot_LC = LC.back()(1);
				if (sdot_cur > sdot_LC)
				{
					result = 1;
					break;
				}
				LC.pop_back();
				tmp_back_list.push_front(Vector2(s_cur, sdot_cur));

				while (true)
				{
					alpha_cur = _topp->determineAlphaBeta(s_cur, sdot_cur)(0);
					_topp->backwardIntegrate(s_cur, sdot_cur, alpha_cur);
					sdot_LC = LC.back()(1);
					if (sdot_cur > sdot_LC)
					{
						result = 1;
						break;
					}
					LC.pop_back();
					tmp_back_list.push_front(Vector2(s_cur, sdot_cur));
					if (LC.size() == 0)
					{
						result = 0;
						break;
					}
				}
				break;
			}
			tmp_back_list.push_front(Vector2(s_cur, sdot_cur));
		}
		s_cur = tmp_back_list.back()(0);
		sdot_cur = tmp_back_list.back()(1);
		LC.insert(LC.end(), tmp_back_list.begin(), tmp_back_list.end());
		return result;
	}

	unsigned int AVP_RRT::backwardIntegrate(Real& s_cur, Real& sdot_cur, std::list<Vector2, Eigen::aligned_allocator<Vector2>>& LC,
		const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints,
		const std::vector<unsigned int>& allMVCPointsFlag)
	{
		unsigned int swi_backward;
		unsigned int swi_backwardVel = 2;
		unsigned int swi_backwardfor = 1;
		while (true)
		{
			if (swi_backwardVel == 2)
				swi_backward = backwardInt(s_cur, sdot_cur, LC, allMVCPoints, allMVCPointsFlag);

			if (swi_backward == 0)
				return false;
			else if (swi_backward == 1)
				break;
			else if (swi_backward == 2)
				swi_backwardVel = backwardIntVel(s_cur, sdot_cur, LC, allMVCPoints, allMVCPointsFlag);

			if (swi_backwardVel == 1)
				break;
			else if (swi_backwardVel == 3)
			{
				swi_backwardfor = backwardforInt(s_cur, sdot_cur, LC, allMVCPoints, allMVCPointsFlag);
				swi_backward = 2;
			}
		}
		return swi_backwardfor;
	}

	unsigned int AVP_RRT::backwardInt(Real& s_cur, Real& sdot_cur, std::list<Vector2, Eigen::aligned_allocator<Vector2>>& LC,
		const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints,
		const std::vector<unsigned int>& allMVCPointsFlag)
	{
		Real alpha_cur, sdot_MVC;
		unsigned int idx;
		int flag;
		while (s_cur > _si)
		{
			alpha_cur = _topp->determineAlphaBeta(s_cur, sdot_cur)(0);
			_topp->backwardIntegrate(s_cur, sdot_cur, alpha_cur);
			if (sdot_cur < 1e-4)
				return 0;

			idx = (unsigned int)round(s_cur / _ds);
			sdot_MVC = allMVCPoints[idx](1);
			flag = allMVCPointsFlag[idx];

			if (sdot_cur > sdot_MVC)
			{
				sdot_cur = sdot_MVC;
				LC.push_front(Vector2(s_cur, sdot_cur));
				if (flag == 2)
					return 1;
				if (flag == 1)
					return 2;
			}
			LC.push_front(Vector2(s_cur, sdot_cur));
		}
		return 1;
	}

	unsigned int AVP_RRT::backwardIntVel(Real& s_cur, Real& sdot_cur, std::list<Vector2, Eigen::aligned_allocator<Vector2>>& LC,
		const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints,
		const std::vector<unsigned int>& allMVCPointsFlag)
	{
		Vector2 alphabeta;
		Real s_next, sdot_next, alpha_cur, beta_cur, slope;
		unsigned int idx;
		int flag;

		while (true)
		{
			s_next = s_cur - _ds;
			if (s_next < _si)
			{
				s_cur = _si; sdot_cur = allMVCPoints[(unsigned int)round(_si / _ds)](1);
				LC.push_front(Vector2(s_cur, sdot_cur));
				return 1;
			}

			idx = (unsigned int)round(s_next / _ds);
			sdot_next = allMVCPoints[idx](1);
			flag = allMVCPointsFlag[idx];
			if (flag == 2)
				return 1;

			alphabeta = _topp->determineAlphaBeta(s_cur, sdot_cur);
			alpha_cur = alphabeta(0);
			beta_cur = alphabeta(1);

			slope = (sdot_next - sdot_cur) / (s_next - s_cur);

			if (-beta_cur <= -slope && -slope <= -alpha_cur)
			{
				s_cur = s_next; sdot_cur = sdot_next;
				LC.push_front(Vector2(s_cur, sdot_cur));
			}
			else if (-alpha_cur < -slope)
				return 2;
			else if (-slope < -beta_cur)
				return 3;
		}
	}

	unsigned int AVP_RRT::backwardforInt(Real& s_cur, Real& sdot_cur, std::list<Vector2, Eigen::aligned_allocator<Vector2>>& LC,
		const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints,
		const std::vector<unsigned int>& allMVCPointsFlag)
	{
		Vector2 alphabeta;
		Real s_next, sdot_next, alpha_cur, beta_cur, slope;
		unsigned int idx;
		int flag;

		while (true)
		{
			s_next = s_cur - _ds;
			if (s_next < _si)
			{
				s_cur = _si; sdot_cur = allMVCPoints[(unsigned int)round(_si / _ds)](1);
				break;
			}
			idx = (unsigned int)round(s_next / _ds);
			sdot_next = allMVCPoints[idx](1);
			flag = allMVCPointsFlag[idx];
			if (flag == 2)
				break;
			alphabeta = _topp->determineAlphaBeta(s_cur, sdot_cur);
			alpha_cur = alphabeta(0);
			beta_cur = alphabeta(1);
			slope = (sdot_next - sdot_cur) / (s_next - s_cur);
			if (-beta_cur <= -slope && -slope <= -alpha_cur)
				break;
			s_cur = s_next; sdot_cur = sdot_next;
		}

		// forward integration
		unsigned int result;
		Real sdot_LC;
		Real s_LC_front = LC.front()(0);
		std::list<Vector2, Eigen::aligned_allocator<Vector2>> tmp_front_list;
		tmp_front_list.push_back(Vector2(s_cur, sdot_cur));
		while (true)
		{
			beta_cur = _topp->determineAlphaBeta(s_cur, sdot_cur)(1);
			_topp->forwardIntegrate(s_cur, sdot_cur, beta_cur);

			if (s_LC_front <= s_cur)
			{
				sdot_LC = LC.front()(1);
				if (sdot_cur > sdot_LC)
				{
					result = 1;
					break;
				}
				LC.pop_front();
				tmp_front_list.push_back(Vector2(s_cur, sdot_cur));

				while (true)
				{
					beta_cur = _topp->determineAlphaBeta(s_cur, sdot_cur)(1);
					_topp->forwardIntegrate(s_cur, sdot_cur, beta_cur);
					sdot_LC = LC.front()(1);
					if (sdot_cur > sdot_LC)
					{
						result = 1;
						break;
					}
					LC.pop_front();
					tmp_front_list.push_back(Vector2(s_cur, sdot_cur));
					if (LC.size() == 0)
					{
						result = 0;
						break;
					}
				}
				break;
			}
			tmp_front_list.push_back(Vector2(s_cur, sdot_cur));
		}
		s_cur = tmp_front_list.front()(0);
		sdot_cur = tmp_front_list.front()(1);
		LC.insert(LC.begin(), tmp_front_list.begin(), tmp_front_list.end());
		return result;
	}

	void AVP_RRT::calculateCLC(std::vector<std::list<Vector2, Eigen::aligned_allocator<Vector2>>>& LC,
		std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& CLC)
	{
		// 0.199 에서 이상해.....
		unsigned int LCsize = LC.size();

		std::cout << "LCsize : " << LCsize << std::endl;

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

		Real tmp = -std::numeric_limits<Real>::max();
		for (unsigned int i = 0; i < LCsize; i++)
		{
			if ((*(it_end[i]))(0) > tmp)
				tmp = (*(it_end[i]))(0);
		}
		s_end = tmp;

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
		std::list<unsigned int> idx_erase;
		Real ds = _topp->getStepSize();
		Real tmp_min;
		unsigned int list_size;

		while (s_cur < s_end)
		{
			idx.clear();

			for (unsigned int i = 0; i < LCsize; i++)
			{
				if (s_cur > ((*(it_end[i]))(0)+_ds*0.5))
				{
					idx_erase.push_back(i);
					LCsize -= 1;
				}
			}

			list_size = idx_erase.size();
			for (unsigned int i = 0; i < list_size; i++)
			{
				it_end.erase(it_end.begin() + idx_erase.front());
				it_begin.erase(it_begin.begin() + idx_erase.front());
				LC.erase(LC.begin() + idx_erase.front());
				idx_erase.pop_front();
				for (std::list<unsigned int>::iterator it = idx_erase.begin(); it != idx_erase.end(); it++)
					*it -= 1;
			}

			for (unsigned int i = 0; i < LCsize; i++)
			{
				if (std::pow(((*(it_begin[i]))(0)- s_cur),2) < _ds*_ds*1e-1) 
					idx.push_back(i);
			}

			LOGIF(((idx.size()) != 0), "CLC function error : idx size must be larger than zero.");

			tmp_min = std::numeric_limits<Real>::max();
			for (unsigned int i = 0; i < idx.size(); i++)
			{
				if ((*(it_begin[idx[i]]))(1) < tmp_min)
					tmp_min = (*(it_begin[idx[i]]))(1);
				it_begin[idx[i]]++;
			}
			CLC.push_back(Vector2(s_cur, tmp_min));

			s_cur += ds;
		}

		idx.clear();
		s_cur = s_end;
		
		for (unsigned int i = 0; i < LCsize; i++)
		{
			if (s_cur >((*(it_end[i]))(0) + _ds*0.5))
			{
				idx_erase.push_back(i);
				LCsize -= 1;
			}
		}

		list_size = idx_erase.size();
		for (unsigned int i = 0; i < list_size; i++)
		{
			it_end.erase(it_end.begin() + idx_erase.front());
			it_begin.erase(it_begin.begin() + idx_erase.front());
			LC.erase(LC.begin() + idx_erase.front());
			idx_erase.pop_front();
			for (std::list<unsigned int>::iterator it = idx_erase.begin(); it != idx_erase.end(); it++)
				*it -= 1;
		}

		for (unsigned int i = 0; i < LCsize; i++)
		{
			if (std::pow(((*(it_begin[i]))(0) - s_cur), 2) < _ds*_ds*1e-1)
				idx.push_back(i);
		}

		tmp_min = std::numeric_limits<Real>::max();
		for (unsigned int i = 0; i < idx.size(); i++)
		{
			if ((*(it_begin[idx[i]]))(1) < tmp_min)
				tmp_min = (*(it_begin[idx[i]]))(1);
			it_begin[idx[i]]++;
		}
		CLC.push_back(Vector2(s_cur, tmp_min));
	}

	unsigned int AVP_RRT::determineAresult(const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints,
		const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& CLC, Real& sdot_beg_star, AVPFLAG avpflag)
	{
		if (CLC.size() == 0)
		{
			if (avpflag == FORWARD)
				sdot_beg_star = allMVCPoints.front()(1);
			else
				sdot_beg_star = allMVCPoints.back()(1);
			return 2;
		}

		bool si_col = false, sf_col = false;

		Real eps_init = std::pow((CLC.front()(0) - _topp->getInitialParam()), 2);
		Real eps_final = std::pow((CLC.back()(0) - _topp->getFinalParam()), 2);

		if (eps_init <= 1e-5)
			si_col = true;
		if (eps_final <= 1e-5)
			sf_col = true;

		if (si_col == false && sf_col == false)
		{
			if (avpflag == FORWARD)
				sdot_beg_star = allMVCPoints.front()(1);
			else // BACKWARD
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
			else // BACKWARD
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
			else // BACKWARD
			{
				sdot_beg_star = CLC.back()(1);
				return 3;
			}
		}
		else// both si_col and sf_col are true
		{
			if (avpflag == FORWARD)
				sdot_beg_star = CLC.front()(1);
			else // BACKWARD
				sdot_beg_star = CLC.back()(1);
			return 5;
		}
	}

	unsigned int AVP_RRT::determineAVPBresult(const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints,
		const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& CLC, const Real sdot_init,
		std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& phi)
	{
		Real s_cur = _si, sdot_cur = sdot_init, beta_cur, sdot_MVC, sdot_CLC;
		bool CLC_Swi = true, CLC_active_swi = false;
		unsigned int cnt = 0;
		int idx;
		
		phi.push_back(Vector2(s_cur, sdot_cur));

		if (CLC.size() == 0)
		{
			CLC_Swi = false;
			sdot_CLC = std::numeric_limits<Real>::max();
		}
			
		while (true)
		{
			beta_cur = _topp->determineAlphaBeta(s_cur, sdot_cur)(1);
			_topp->forwardIntegrate(s_cur, sdot_cur, beta_cur);
			idx = (unsigned int)round(s_cur / _ds);
			sdot_MVC = allMVCPoints[idx](1);
			
			if (CLC_Swi)
				if ((CLC.front()(0) - s_cur) < _ds*0.1)
					CLC_active_swi = true;

			if (CLC_active_swi)
			{
				if (cnt < CLC.size())
					sdot_CLC = CLC[cnt++](1);
				else
					sdot_CLC = std::numeric_limits<Real>::max();
			}
				
			if (sdot_cur < 0)
				return 1;
			else if (s_cur > _sf)
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
		int idx;
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
			idx = (unsigned int)round(s_cur / ds);
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
		int idx;
		si = _topp->getInitialParam();
		ds = _topp->getStepSize();
		s_cur = _topp->getFinalParam();
		sdot_cur = sdot_test;

		int phiSize = phi.size();

		while (true)
		{
			alpha_cur = _topp->determineAlphaBeta(s_cur, sdot_cur)(0);
			_topp->backwardIntegrate(s_cur, sdot_cur, alpha_cur);
			idx = (unsigned int)round(s_cur / ds);
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
		int idx;
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
			idx = (unsigned int)round(s_cur / ds);
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
		bool isValid = true;

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