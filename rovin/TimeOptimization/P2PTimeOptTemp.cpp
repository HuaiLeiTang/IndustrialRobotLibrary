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

	void AVP_RRT::generateTrajectorySegment(int idx)
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
			extendTree(&_startTree, qrand);
			extendTree(&_goalTree, qrand);

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
		//make random config and save at qrand;
	}

	Vertex * AVP_RRT::extendTree(Tree * tree, const VectorX qrand) // 반환 Vertex* 로??
	{
		double tmpDist;
		Vertex * nVertex = tree->findNearestNeighbor(qrand, &tmpDist);

		Vector2 endInterval;
		Vertex * candVertex = runAVP(nVertex, endInterval);


		return candVertex;

	}

	Vertex * AVP_RRT::runAVP(Vertex * nVertex, Vector2 & endInterval)
	{
		///// nVertex : nearest vertex, 아래 주석처럼 사용.
		//Vector2 begInterval = nVertex->_interval;
		//MatrixX inpath = nVertex->_inpath;

		endInterval(0) = 0.0;
		endInterval(1) = 0.5;

		Vertex * candVertex = new Vertex();

		// if succeeded avp.
		// candVertex setting 함수들 불러서 멤버변수4개 세팅해주고 return;
		return candVertex;
	}




	void Tree::initializeTree(Vertex * rootVertex) {
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