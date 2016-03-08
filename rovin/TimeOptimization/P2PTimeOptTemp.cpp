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
		// start tree, goal tree -> waypont 로 initialize
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

			// 이렇게 하면 안되지롱... avp_backward 따로 있으니까... 나눠서.... 흠 일단 ㄱ
			// extendTree도 vertex 받기.....
			Vertex * candiVertex = extendTree(&_startTree, qrand, true);
			if (candiVertex != NULL)
			{
				// add vertex
				_startTree.addVertex(candiVertex);
				
				// connection test to the _goalTree
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
			startVertex = new Vertex(_waypoints[idx].getJointq(), Vector2().setZero(), MatrixX().setZero(), NULL);
		else
			startVertex = new Vertex(_waypoints[idx].getJointq(), _wayPointInterval, MatrixX().setZero(), NULL);

		Vertex * goalVertex = new Vertex(_waypoints[idx + 1].getJointq(), Vector2().setZero(), MatrixX().setZero(), NULL);


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
		MatrixX Pnew;
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
		for (int i = 0; i < Pnew.cols(); i++)
		{
			// TO DO:
			// collision check! of Pnew.col(i)
			// if coliision detected, isSucceeded = false and break!
		}

		if (isSucceeded)
			return new Vertex(qnew, endInterval, Pnew, nVertex);
		else
			return NULL;

	}

	void AVP_RRT::interpolate(Vertex * nVertex, const VectorX qrand, const double dist, MatrixX & Pnew, VectorX & qnew)
	{
		if (dist < _stepsize)
			qnew = qrand;
		else
			qnew = nVertex->_config + (qrand - nVertex->_config)*_stepsize / dist;

		// TO DO: INTERPOLATE


	}

	bool AVP_RRT::testConnection(Vertex * vertex, Tree * tree, Vertex ** oVertex)
	{
		// tree 에서 vertex가 extended 되었을때, vertex랑 tree를 집어넣어서 connection test..
		// oVertex 는 tree에서 연결이 되면 vertex point 아니면... NULL로 가자
		// 연결되면 treu, 안되면 false..
		// 이 함수에서도 avp test 해야함

		return false;
	}

	bool AVP_RRT::runAVP(const MatrixX& Pnew, const Vector2& nearInterval, Vector2 & endInterval)
	{
		endInterval(0) = 0;
		endInterval(1) = 0;
		return false;
	}

	bool AVP_RRT::runAVPbackward(const MatrixX& Pnew, const Vector2& nearInterval, Vector2 & endInterval)
	{
		endInterval(0) = 0;
		endInterval(1) = 0;
		return false;
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