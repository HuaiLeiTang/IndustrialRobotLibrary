#pragma once

#include "TOPP.h"

// inpath �� bspline ���� ? �ƴϸ� discrete �ϰ�??? �׸��� inpath�� q0 ~ q1 ������ �����ΰ� �׸��� linear �ϰ� ����� �ǳ�?
// s, sdot �� �������ʿ����??
// RRT step size ������ ��� ������ �־���ϳ���?????


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
		AVP_RRT(const SerialOpenChainPtr& robot);

		void setWayPoints(const std::vector<WayPoint>& waypoints) { _waypoints = waypoints; }
		void addWayPoints(const WayPoint& waypoint) { _waypoints.push_back(waypoint); }
		void clearWayPoints() { _waypoints.clear(); }

		const std::vector<WayPoint>& getWayPoints() const { return _waypoints; }
		const std::vector<Tree*>& getTree() const { return _tree; }

		void generateTrajectory();

	private:
		void initialization();
		void makeRandomConfig(VectorX& qrand);
		

	private:
		SerialOpenChainPtr _robot;
		
		TOPP* topp;

		std::vector<WayPoint> _waypoints;
		std::vector<Tree*> _tree;


		BSpline<-1, -1, -1> _finalPath;
		// MatrixX _finalPath; //---> �̰� �� �����ƺ���

		unsigned int _dof;
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
		BSpline<-1, -1, -1> _inpath; // VectorX _inpath; ---> ������.... spline ��ߵǳ�??
		// MatrixX _inpath; // --> �̰� �� ���ƺ���
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

		bool extend(VectorX& qrand);
		void makepath(); // ���� path ����

	private:
		AVP_RRT* _avp_rrt;
		Extend* _extend;
		std::vector<Vertex> _vertices;
		
		BSpline<-1, -1, -1> _path; // VectorX _inpath; ---> ������.... spline ��ߵǳ�??
		// MatrixX _path; // --> �̰� �� ���ƺ���

		WayPoint _initPoint;
		WayPoint _finalPoint;
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

	private:


	private:
		Tree* _tree;
		
		Vertex _Vnear;
		Vertex _Vnew;

		VectorX _qrand;
		VectorX _qnew;
		Vector2 _intervalnew;
		BSpline<-1, -1, -1> _pathnew; // VectorX _inpath; ---> ������.... spline ��ߵǳ�??
		unsigned int _nearestVertexIdx;
	};
}