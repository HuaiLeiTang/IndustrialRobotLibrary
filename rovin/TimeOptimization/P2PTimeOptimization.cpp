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
		delete topp;
	}
	
	AVP_RRT::AVP_RRT(const SerialOpenChainPtr & robot, CONSTRAINT_TYPE constraintType)
	{
		_robot = robot;
		_dof = robot->getNumOfJoint();

		Real ds = 1e-3, vi = 0, vf = 0, si = 0, sf = 1;
		topp = new TOPP(_robot, vi, vf, ds, si, sf, constraintType);
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

		// collisioncheck �ؼ� collision �Ǹ� ��� �ٽ� ����? --> path�� �ٽ� ������ �ǳ���??
		// AVP ���� timeparameteriable ���� ������ ��� �ٽ� ����?
	}

	void Tree::makepath()
	{
	}

	/////////////////// Extend class ///////////////////

	void Extend::update()
	{
		_Vnew.setconfig(_qnew);
		//_Vnew.setinpath(_pathnew);
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
		
		// make topp
		_tree->_avp_rrt->topp->setJointTrajectory(_pathnew);

		// Step A. Computing the limiting curves
		// Acalculate MVC
		std::vector<Vector2> allMVCPoints;
		_tree->_avp_rrt->topp->calculateAllMVCPoint();
		allMVCPoints = _tree->_avp_rrt->topp->getAllMVCPoint();

		// calculate CLC
		std::vector<SwitchPoint> allSwitchPoint;
		_tree->_avp_rrt->topp->calculateAllSwitchPoint();
		allSwitchPoint = _tree->_avp_rrt->topp->getAllSwitchPoint();

		// A1


		// A2-5


		return true;
	}

	

}