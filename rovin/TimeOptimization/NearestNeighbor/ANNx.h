#ifndef ANNx_H
#define ANNx_H

#include "ANN.h"			// ANN includes

enum	{ANN_LO=0, ANN_HI=1};	// splitting indices
enum	{ANN_IN=0, ANN_OUT=1};	// shrinking indices

enum	{LO=0, HI=1};			// splitting indices
enum	{IN=0, OUT=1};			// shrinking indices

// what to do in case of error
enum ANNerr {ANNwarn = 0, ANNabort = 1};

extern int		ANNmaxPtsVisited;	// maximum number of pts visited
extern int		ANNptsVisited;		// number of pts visited in search

class ANNorthRect {
public:
	ANNpoint		lo;			// rectangle lower bounds
	ANNpoint		hi;			// rectangle upper bounds
	//
	ANNorthRect(				// basic constructor
		int				dd,			// dimension of space
		ANNcoord		l=0,		// default is empty
		ANNcoord		h=0)
	{  lo = annAllocPt(dd, l);  hi = annAllocPt(dd, h); }

	ANNorthRect(				// (almost a) copy constructor
		int				dd,			// dimension
		const			ANNorthRect &r) // rectangle to copy
	{  lo = annCopyPt(dd, r.lo);  hi = annCopyPt(dd, r.hi);  }

	ANNorthRect(				// construct from points
		int				dd,			// dimension
		ANNpoint		l,			// low point
		ANNpoint		h)			// hight point
	{  lo = annCopyPt(dd, l);  hi = annCopyPt(dd, h);  }

	~ANNorthRect()				// destructor
	{  annDeallocPt(lo);  annDeallocPt(hi);  }

	ANNbool inside(int dim, ANNpoint p);// is point p inside rectangle?
};

void annAssignRect(				// assign one rect to another
				   int				dim,		// dimension (both must be same)
				   ANNorthRect		&dest,		// destination (modified)
				   const ANNorthRect &source);	// source

//----------------------------------------------------------------------
//	Orthogonal (axis aligned) halfspace
//	An orthogonal halfspace is represented by an integer cutting
//	dimension cd, coordinate cutting value, cv, and side, sd, which is
//	either +1 or -1. Our convention is that point q lies in the (closed)
//	halfspace if (q[cd] - cv)*sd >= 0.
//----------------------------------------------------------------------

class ANNorthHalfSpace {
public:
	int				cd;			// cutting dimension
	ANNcoord		cv;			// cutting value
	int				sd;			// which side
	//
	ANNorthHalfSpace()			// default constructor
	{  cd = 0; cv = 0;  sd = 0;  }

	ANNorthHalfSpace(			// basic constructor
		int				cdd,		// dimension of space
		ANNcoord		cvv,		// cutting value
		int				sdd)		// side
	{  cd = cdd;  cv = cvv;  sd = sdd;  }

	ANNbool in(ANNpoint q) const	// is q inside halfspace?
	{  return  (ANNbool) ((q[cd] - cv)*sd >= 0);  }

	ANNbool out(ANNpoint q) const	// is q outside halfspace?
	{  return  (ANNbool) ((q[cd] - cv)*sd < 0);  }

	ANNdist dist(ANNpoint q) const	// (squared) distance from q
	{  return  (ANNdist) ANN_POW(q[cd] - cv);  }

	void setLowerBound(int d, ANNpoint p)// set to lower bound at p[i]
	{  cd = d;  cv = p[d];  sd = +1;  }

	void setUpperBound(int d, ANNpoint p)// set to upper bound at p[i]
	{  cd = d;  cv = p[d];  sd = -1;  }

	void project(ANNpoint &q)		// project q (modified) onto halfspace
	{  if (out(q)) q[cd] = cv;  }
};

// array of halfspaces
typedef ANNorthHalfSpace *ANNorthHSArray;

#endif
