#pragma once

#include "TOPP.h"

// inpath 는 bspline 으로 ? 아니면 discrete 하게??? 그리고 inpath는 q0 ~ q1 사이의 값들인가 그리고 linear 하게 만들면 되나?
// RRT step size 같은건 어디서 가지고 있어야하나용?????
// TOPP 는 누가 가지고 있는게 좋을까 ----> 가장 큰 클래스!

namespace rovin
{
	class AVP_RRT;
	class Vertex;
	class Tree;
	class Extend;
	class WayPoint;

	class AVP_RRT
	{
		friend class Tree;
		friend class Extend;
	public:
		AVP_RRT();
		~AVP_RRT();
		AVP_RRT(const SerialOpenChainPtr& robot, CONSTRAINT_TYPE constraintType = CONSTRAINT_TYPE::TORQUE);

		void setWayPoints(const std::vector<WayPoint>& waypoints) { _waypoints = waypoints; }
		void addWayPoints(const WayPoint& waypoint) { _waypoints.push_back(waypoint); }
		void clearWayPoints() { _waypoints.clear(); }

		const std::vector<WayPoint>& getWayPoints() const { return _waypoints; }
		const std::vector<Tree*>& getTree() const { return _tree; }
		void generateTrajectory();

	public:
		void initialization();
		void makeRandomConfig(VectorX& qrand);
		

	public:
		


		SerialOpenChainPtr _robot;
		TOPP* _topp;
		std::vector<WayPoint> _waypoints;
		std::vector<Tree*> _tree;

		BSpline<-1, -1, -1> _finalPath;
		// MatrixX _finalPath; //---> 이게 더 괜찮아보임
		unsigned int _dof;



	public:
		//////// function for test ///////
		void setTOPP(TOPP* topp) { _topp = topp; }
		void extendsetting() {}
		void addtree(Tree* tree) { _tree.push_back(tree); }
	};

	class Vertex
	{
	public:
		Vertex() : _config(VectorX()), _interval(Vector2()), _inpath(BSpline<-1, -1, -1>()), _parentVertexIdx(-1) {}
		~Vertex() {}
		Vertex(const VectorX& config, const Vector2& interval, const BSpline<-1, -1, -1>& inpath, const int parentVertexIdx)
			: _config(config), _interval(interval), _inpath(inpath), _parentVertexIdx(parentVertexIdx) {}

		void setconfig(const VectorX& config) { _config = config; }
		void setinterval(const Vector2& interval) { _interval = interval; }
		void setinpath(const BSpline<-1, -1, -1>& inpath) { _inpath = inpath; }
		void setparentVertexIdx(const int parentVertexIdx) { _parentVertexIdx = parentVertexIdx; }

	private:
		VectorX _config;
		Vector2 _interval;
		BSpline<-1, -1, -1> _inpath; // VectorX _inpath; ---> 생각좀.... spline 써야되나??
		// MatrixX _inpath; // --> 이게 더 좋아보임
		int _parentVertexIdx;
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

	class Tree
	{
		friend class AVP_RRT;
		friend class Extend;
	public:
		Tree() {}
		~Tree() { if (_extend != NULL) delete _extend; }
		Tree(AVP_RRT* avp_rrt, const WayPoint& initPoint, const WayPoint& finalPoint);
		Tree(AVP_RRT* avp_rrt) : _avp_rrt(avp_rrt){}

		bool extend(VectorX& qrand);
		void makepath(); // 최종 path 생성

	public:
		AVP_RRT* _avp_rrt;
		Extend* _extend;
		std::vector<Vertex> _vertices;
		
		BSpline<-1, -1, -1> _path; // VectorX _inpath; ---> 생각좀.... spline 써야되나??
		// MatrixX _path; // --> 이게 더 좋아보임

		WayPoint _initPoint;
		WayPoint _finalPoint;

	public:
		//////// function for test ///////
		void addextend(Extend* extend) { _extend = extend; }

	};

	class Extend
	{
		friend class AVP_RRT;
		friend class Tree;
	public:
		Extend() {}
		~Extend() {}
		Extend(Tree* tree) : _tree(tree) {}
		Extend(Tree* tree, const VectorX& qrand) : _tree(tree), _qrand(qrand) {}

		void update();
		void nearestNeighbor(VectorX& qrand);
		void interpolate();
		bool collisioncheck();
		bool AVP();

	public:
		bool calculateLimitingCurves(const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints, const std::vector<SwitchPoint>& allSwitchPoint, std::vector<std::list<Real>>& LC);

	public:
		Tree* _tree;
		
		Vertex _Vnear;
		Vertex _Vnew;

		VectorX _qrand;
		VectorX _qnew;
		Vector2 _intervalnew;
		MatrixX _pathnew; //BSpline<-1, -1, -1> _pathnew; // VectorX _inpath; ---> 생각좀.... spline 써야되나??
		unsigned int _nearestVertexIdx;
	};
}