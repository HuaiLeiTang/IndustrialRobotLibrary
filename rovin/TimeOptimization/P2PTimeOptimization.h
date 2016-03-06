#pragma once

#include "TOPP.h"


namespace rovin
{
	class AVP_RRT;
	class Extend;

	class Tree;
	class Vertex;
	class WayPoint;



	class AVP_RRT
	{
	public:
		AVP_RRT() {}
		~AVP_RRT() {}
		AVP_RRT(const SerialOpenChainPtr& robot);

		void setWayPoints(const std::vector<WayPoint>& waypoints) { _waypoints = waypoints; }
		void addWayPoints(const WayPoint& waypoint) { _waypoints.push_back(waypoint); }
		void clearWayPoints() { _waypoints.clear(); }

		const std::vector<WayPoint>& getWayPoints() const { return _waypoints; }

	private:
		void makeRandomConfig(VectorX& qrand);


	private:
		SerialOpenChainPtr _robot;
		Extend* _extend;
		TOPP* topp;

		std::vector<WayPoint> _waypoints;
		std::vector<Tree> _tree;


		unsigned int _dof;
	};


	class Tree
	{
	public:
		Tree() {}
		~Tree() {}

	private:
		// Vertex 를 가지고 있어야..

	};

	class Vertex
	{
	public:
		Vertex() {}
		~Vertex() {}

	private:
		VectorX _config;
		// VectorX _inpath;
		Vector2 _interval;


	};

	class WayPoint
	{
	public:
		WayPoint() : _q(VectorX()), _qdot(VectorX()) {}
		~WayPoint() {}
	private:
		void setJointq(const VectorX& q) { _q = q; }
		void setJointqdot(const VectorX& qdot) { _qdot = qdot; }

		const VectorX& getJointq() const { return _q; }
		const VectorX& getJointqdot() const { return _qdot; }
	private:
		VectorX _q;
		VectorX _qdot;
	};

}