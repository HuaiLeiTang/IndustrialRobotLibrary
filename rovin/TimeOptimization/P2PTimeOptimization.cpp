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
	}
	
	AVP_RRT::AVP_RRT(const SerialOpenChainPtr & robot)
	{
		_robot = robot;
		_dof = robot->getNumOfJoint();
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
		_Vnew.setinpath(_pathnew);
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

		return true;
	}

	

}