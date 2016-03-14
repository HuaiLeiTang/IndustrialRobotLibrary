#pragma once

#include "TOPP.h"
#include "NearestNeighbor/multiann.h"

#include <fstream>
#include <string>

// inpath 는 bspline 으로 ? 아니면 discrete 하게??? 그리고 inpath는 q0 ~ q1 사이의 값들인가 그리고 linear 하게 만들면 되나?
// s, sdot 은 저장할필요없나?? 없을듯
// RRT step size 같은건 어디서 가지고 있어야하나용?????
// TOPP 는 누가 가지고 있는게 좋을까 ----> 가장 큰 클래스!!

#define MPNN_TREE_MAXIMUM_NODES	100000
#define AVP_RRT_MAX_ITER		100000

namespace rovin
{
	class Vertex;
	class Tree;
	class AVP_RRT;
	class WayPoint;

	class Vertex
	{
		friend class Tree;
		friend class AVP_RRT;
	public:
		Vertex() : _config(VectorX()), _interval(Vector2()), _inpath(std::list<VectorX>()), _parentVertex(NULL) {}
		~Vertex() { }
		Vertex(const VectorX& config, const Vector2& interval, const std::list<VectorX>& inpath, Vertex * parentVertex)
			: _config(config), _interval(interval), _inpath(inpath), _parentVertex(parentVertex) {}

		void setconfig(const VectorX& config) { _config = config; }
		void setinterval(const Vector2& interval) { _interval = interval; }
		void setinpath(const std::list<VectorX>& inpath) { _inpath = inpath; }
		void setparentVertex(Vertex * parentVertex) { _parentVertex = parentVertex; }

	private:
		VectorX _config;
		Vector2 _interval;
		//MatrixX _inpath;
		std::list<VectorX> _inpath;
		Vertex * _parentVertex;
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	};

	class Tree
	{
		friend class AVP_RRT;
	public:
		// robotdof 대신 rootVertex->_config.size() 해도 될듯..
		Tree() {}
		//Tree(Vertex * rootVertex) {
		//	initializeTree(rootVertex);
		//}
		~Tree() {
			clearTree();
		}

		void clearTree()
		{
			if (_mpnnTree != NULL)
				delete _mpnnTree;
			for (unsigned int i = 0; i < _nodes.size(); i++)
				if (_nodes[i] != NULL)
					delete _nodes[i];
			_nodes.clear();

		}

		void initializeTree(Vertex * rootVertex);

		Vertex * findNearestNeighbor(VectorX targetConfig, double * distance);

		// addVertex 부르기 전에 vertex 멤버변수 다 값 저장한다음에 넣기
		void addVertex(Vertex* vertex);

	private:
		MultiANN * _mpnnTree;
		//Vertex * _rootVertex;
		std::vector<Vertex*> _nodes; ///> [0] component has root vertex;
									 //unsigned int _numNodes;

	};

	class AVP_RRT
	{
		
	public:
		enum RETURNFLAG {SUCCESS = 1, EXCEED_MAX_ITER = -1, EXCEED_MAX_NODE = -2}; // more flags are needed

	public:
		AVP_RRT() {
			srand(time(NULL));
		};
		~AVP_RRT() {}
		AVP_RRT(const SerialOpenChainPtr& robot, CONSTRAINT_TYPE constraintType){
			_robot = robot;
			_constraintType = constraintType; 
			
			Real ds = 1e-3, vi = 0, vf = 0, si = 0, sf = 1;
			_topp = TOPPPtr(new TOPP(_robot, vi, vf, ds, si, sf, constraintType));

			_dof = robot->getNumOfJoint();
			_stepsize = 0.1;
			_curSegment = 0;

			srand(time(NULL));
		}

		void setWayPoints(const std::vector<WayPoint>& waypoints) { _waypoints = waypoints; }
		void addWayPoints(const WayPoint& waypoint) { _waypoints.push_back(waypoint); }
		void clearWayPoints() { _waypoints.clear(); }

		const std::vector<WayPoint>& getWayPoints() const {	return _waypoints; }
		const MatrixX& getFinalPath() const { return _finalPath; }

		// int 받아서 i-th segment에 대해서 도는 generate Trajectory 만들고, 최종 경로를 저장하는 컨테이너 하나 만들기.. vector로?
		// idx는 waypoint 개수 -1
		void generateTrajectory();
		RETURNFLAG generateTrajectorySegment(int idx);

		void treeInitialization(int idx);


	private:
		void makeRandomConfig(VectorX& qrand);
		Vertex * extendTree(Tree* tree, const VectorX qrand, bool atStartTree);
		void interpolate(Vertex * nVertex, const VectorX qrand, const double dist, /* OUTPUT */ std::list<VectorX>& Pnew, VectorX& qnew);
		bool testConnection(Vertex * vertex, Tree * tree, bool forward, /* OUTPUT */ Vertex ** cVertex, Vertex ** oVertex);
		bool checkIntersectionOfTwoIntervals(const Vector2& vec1, const Vector2& vec2);
		void extractPath(Vertex * sVertex, Vertex * gVertex, int idx);

		// calculate for AVP
	public:
		bool runAVP(std::list<VectorX>& Pnew, Vector2& nearInterval, /* OUTPUT */ Vector2& endInterval);
		bool runAVPbackward(std::list<VectorX>& Pnew, Vector2& nearInterval, /* OUTPUT */ Vector2& endInterval);
		bool calculateLimitingCurves(const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints,
			const std::vector<unsigned int>& allMVCPointsFlag, const std::vector<SwitchPoint>& allSwitchPoint,
			std::vector<std::list<Vector2, Eigen::aligned_allocator<Vector2>>>& LC);
		void calculateCLC(std::vector<std::list<Vector2, Eigen::aligned_allocator<Vector2>>>& LC,
			std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& CLC);

	private:
		SerialOpenChainPtr _robot;

		CONSTRAINT_TYPE _constraintType;
		std::vector<WayPoint> _waypoints; // contains q and qdot
		
		Tree _startTree;
		Tree _goalTree; ///> start and goal trees are made for every segment
		std::vector<std::list<VectorX>> _segmentPath;
		MatrixX _finalPath; ///> concatenation of _segmentPath

		TOPPPtr _topp;
		unsigned int _dof;
		unsigned int _numSegment;
		unsigned int _curSegment;
		double _stepsize;
		Vector2 _wayPointInterval; // useless...?


		//////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////

	public:
		std::vector<Real> s;
		std::vector<Real> sdot_MVC;

		std::vector<Real> s_LC;
		std::vector<Real> sdot_LC;
		std::vector<Real> sdot_CLC;

		void saveRealVector2txt(std::vector<Real> in, std::string filename)
		{
			std::ofstream fout;
			fout.open(filename);

			for (unsigned int i = 0; i < in.size(); i++)
				fout << in[i] << std::endl;

			fout.close();
		}

		void saveMVC(std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints)
		{
			std::string name = "C:/Users/crazy/Desktop/Time optimization";
			std::vector<Real> s, sdot;
			for (unsigned int i = 0; i < allMVCPoints.size(); i++)
			{
				s.push_back(allMVCPoints[i](0));
				sdot.push_back(allMVCPoints[i](1));
			}

			saveRealVector2txt(s, "C:/Users/crazy/Desktop/Time optimization/s.txt");
			saveRealVector2txt(sdot, "C:/Users/crazy/Desktop/Time optimization/sdot_MVC.txt");
		}

		void saveRealList2txt(std::list<Real> in, std::string filename)
		{
			std::vector<Real> vec;
			for (std::list<Real>::iterator it = in.begin(); it != in.end(); ++it)
			{
				vec.push_back(*(it));
			}
			saveRealVector2txt(vec, filename);
		}

		void saveLC(std::list<Vector2, Eigen::aligned_allocator<Vector2>>& in, std::string filename_s, std::string filename_sdot)
		{
			std::vector<Real> s;
			std::vector<Real> sdot;

			Vector2 tmp;

			for (std::list<Vector2, Eigen::aligned_allocator<Vector2>>::iterator it = in.begin(); it != in.end(); ++it)
			{
				tmp = *(it);
				s.push_back(tmp[0]);
				sdot.push_back(tmp[1]);
			}
			saveRealVector2txt(s, filename_s);
			saveRealVector2txt(sdot, filename_sdot);
		}


	};

	class WayPoint
	{
	public:
		WayPoint(const VectorX& q, const VectorX& qdot) : _q(q), _qdot(qdot) {}
		WayPoint() {}
		~WayPoint() {}
	public:
		void setJointq(const VectorX& q) { _q = q; }
		void setJointqdot(const VectorX& qdot) { _qdot = qdot; }
		const VectorX& getJointq() const { return _q; }
		const VectorX& getJointqdot() const { return _qdot; }
	private:
		VectorX _q;
		VectorX _qdot;
	};

}