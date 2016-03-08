#ifndef ANN_H
#define ANN_H



#include <cmath>			// math includes
#include <iostream>			// I/O streams

#ifdef ANN_NO_LIMITS_H					// limits.h unavailable
#include <cvalues>					// replacement for limits.h
const double ANN_DBL_MAX = MAXDOUBLE;	// insert maximum double
#else
#include <climits>
#include <cfloat>
const double ANN_DBL_MAX = DBL_MAX;
#endif


enum ANNbool {ANNfalse = 0, ANNtrue = 1}; // ANN boolean type (non ANSI C++)

typedef double	ANNcoord;				// coordinate data type
typedef double	ANNdist;				// distance data type

typedef int		ANNidx;					// point index
const ANNidx	ANN_NULL_IDX = -1;		// a NULL point index

const ANNdist	ANN_DIST_INF = ANN_DBL_MAX;


#ifdef DBL_DIG							// number of sig. bits in ANNcoord
const int	 ANNcoordPrec	= DBL_DIG;
#else
const int	 ANNcoordPrec	= 15;	// default precision
#endif

const ANNbool	ANN_ALLOW_SELF_MATCH	= ANNtrue;

#define ANN_POW(v)			((v)*(v))
#define ANN_ROOT(x)			sqrt(x)
#define ANN_SUM(x,y)		((x) + (y))
#define ANN_DIFF(x,y)		((y) - (x))
#define ANN_MIN(x,y)		(x > y) ? y : x


typedef ANNcoord* ANNpoint;			// a point
typedef ANNpoint* ANNpointArray;	// an array of points 
typedef ANNdist*  ANNdistArray;		// an array of distances 
typedef ANNidx*   ANNidxArray;		// an array of point indices



ANNdist annDist(
				int				dim,		// dimension of space
				ANNpoint		p,			// points
				ANNpoint		q);

ANNpoint annAllocPt(
					int				dim,		// dimension
					ANNcoord		c = 0);		// coordinate value (all equal)

ANNpointArray annAllocPts(
						  int				n,			// number of points
						  int				dim);		// dimension

void annDeallocPt(
				  ANNpoint		&p);		// deallocate 1 point

void annDeallocPts(
				   ANNpointArray	&pa);		// point array

ANNpoint annCopyPt(
				   int				dim,		// dimension
				   ANNpoint		source);	// point to copy

class ANNpointSet {
public:
	virtual ~ANNpointSet() {}			// virtual distructor

	virtual void annkSearch(			// approx k near neighbor search
		ANNpoint		q,				// query point
		int				k,				// number of near neighbors to return
		ANNidxArray		nn_idx,			// nearest neighbor array (modified)
		ANNdistArray	dd,				// dist to near neighbors (modified)
		double			eps=0.0			// error bound
		) = 0;							// pure virtual (defined elsewhere)

	virtual int theDim() = 0;			// return dimension of space
	virtual int nPoints() = 0;			// return number of points
	// return pointer to points
	virtual ANNpointArray thePoints() = 0;
};



class ANNbruteForce: public ANNpointSet {
	int				dim;				// dimension
	int				n_pts;				// number of points
	ANNpointArray	pts;				// point array
public:
	ANNbruteForce(						// constructor from point array
		ANNpointArray	pa,				// point array
		int				n,				// number of points
		int				dd);			// dimension

	~ANNbruteForce();					// destructor

	void annkSearch(					// approx k near neighbor search
		ANNpoint		q,				// query point
		int				k,				// number of near neighbors to return
		ANNidxArray		nn_idx,			// nearest neighbor array (modified)
		ANNdistArray	dd,				// dist to near neighbors (modified)
		double			eps=0.0);		// error bound

	int theDim()						// return dimension of space
	{ return dim; }

	int nPoints()						// return number of points
	{ return n_pts; }

	ANNpointArray thePoints()			// return pointer to points
	{  return pts;  }
};


enum ANNsplitRule {
	ANN_KD_STD				= 0,	// the optimized kd-splitting rule
	ANN_KD_MIDPT			= 1,	// midpoint split
	ANN_KD_FAIR				= 2,	// fair split
	ANN_KD_SL_MIDPT			= 3,	// sliding midpoint splitting method
	ANN_KD_SL_FAIR			= 4,	// sliding fair split method
	ANN_KD_SUGGEST			= 5};	// the authors' suggestion for best
	const int ANN_N_SPLIT_RULES		= 6;	// number of split rules


	class ANNkd_node;				// generic node in a kd-tree
	typedef ANNkd_node*	ANNkd_ptr;	// pointer to a kd-tree node

	class ANNkd_tree: public ANNpointSet {
	protected:
		int				dim;				// dimension of space
		int				n_pts;				// number of points in tree
		int				bkt_size;			// bucket size
		ANNpointArray	pts;				// the points
		ANNidxArray		pidx;				// point indices (to pts array)
		ANNkd_ptr		root;				// root of kd-tree
		ANNpoint		bnd_box_lo;			// bounding box low point
		ANNpoint		bnd_box_hi;			// bounding box high point

		double			*Scale;				// scaling
		int				*Topology;			// topology

		void SkeletonTree(					// construct skeleton tree
			int				n,				// number of points
			int				dd,				// dimension
			int				bs,				// bucket size
			ANNpointArray pa = NULL,		// point array (optional)
			ANNidxArray pi = NULL);			// point indices (optional)

	public:
		ANNkd_tree(							// build skeleton tree
			int				n = 0,			// number of points
			int				dd = 0,			// dimension
			int				bs = 1);		// bucket size

		ANNkd_tree(							// build from point array
			ANNpointArray	pa,				// point array
			int				n,				// number of points
			int				dd,				// dimension
			double	       	*ss,			// scaling 
			int			 	*tt,			// topology of space 
			int				bs = 1,			// bucket size
			ANNsplitRule	split = ANN_KD_SUGGEST);	// splitting method

		~ANNkd_tree();						// tree destructor

		void annkSearch(					// approx k near neighbor search
			ANNpoint		q,				// query point
			int				k,				// number of near neighbors to return
			ANNidxArray		nn_idx,			// nearest neighbor array (modified)
			ANNdistArray	dd,				// dist to near neighbors (modified)
			double			eps=0.0);		// error bound

		int theDim()						// return dimension of space
		{ return dim; }

		int nPoints()						// return number of points
		{ return n_pts; }

		ANNpointArray thePoints()			// return pointer to points
		{  return pts;  }

	};								

	void annMaxPtsVisit(	// max. pts to visit in search
		int				maxPts);	// the limit

	void annClose();		// called to end use of ANN

#endif
