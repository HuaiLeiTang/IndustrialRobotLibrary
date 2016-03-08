#ifndef ANN_kd_util_H
#define ANN_kd_util_H

#include "kd_tree.h"					// kd-tree declarations

//----------------------------------------------------------------------
//	externally accessible functions
//----------------------------------------------------------------------

double annAspectRatio(			// compute aspect ratio of box
					  int					dim,			// dimension
					  const ANNorthRect	&bnd_box);		// bounding cube

void annEnclRect(				// compute smallest enclosing rectangle
				 ANNpointArray		pa,				// point array
				 ANNidxArray			pidx,			// point indices
				 int					n,				// number of points
				 int					dim,			// dimension
				 ANNorthRect &bnds);					// bounding cube (returned)

void annEnclCube(				// compute smallest enclosing cube
				 ANNpointArray		pa,				// point array
				 ANNidxArray			pidx,			// point indices
				 int					n,				// number of points
				 int					dim,			// dimension
				 ANNorthRect &bnds);					// bounding cube (returned)

class ANNRectangle
{
	ANNcoord     	bnds1[2];		// bounding values
	ANNcoord      bnds2[2];		// of the rectangle
	ANNcoord     	bnds3[2];		// in four dimensions
	ANNcoord     	bnds4[2];
public:
	ANNRectangle(			// constructor
		ANNcoord lb1, ANNcoord hb1,
		ANNcoord lb2, ANNcoord hb2,
		ANNcoord lb3, ANNcoord hb3,
		ANNcoord lb4, ANNcoord hb4 );

	~ANNRectangle();			// destructor

	ANNdist DistPointRectangle (ANNcoord p1, ANNcoord p2, ANNcoord p3, ANNcoord p4);      
};


ANNdist annBoxDistance(			// compute distance from point to box
					   const ANNpoint		q,				// the point
					   const ANNpoint		lo,				// low point of box
					   const ANNpoint		hi,				// high point of box
					   int					dim,			// dimension of space	
					   double 				*Scale,			// scaling coefficients
					   int 				*Topology);		// topology of space

ANNcoord annSpread(				// compute point spread along dimension
				   ANNpointArray		pa,				// point array
				   ANNidxArray			pidx,			// point indices
				   int					n,				// number of points
				   int					d);				// dimension to check

void annMinMax(					// compute min and max coordinates along dim
			   ANNpointArray		pa,				// point array
			   ANNidxArray			pidx,			// point indices
			   int					n,				// number of points
			   int					d,				// dimension to check
			   ANNcoord&			min,			// minimum value (returned)
			   ANNcoord&			max);			// maximum value (returned)

int annMaxSpread(				// compute dimension of max spread
				 ANNpointArray		pa,				// point array
				 ANNidxArray			pidx,			// point indices
				 int					n,				// number of points
				 int					dim);			// dimension of space

void annMedianSplit(			// split points along median value
					ANNpointArray		pa,				// points to split
					ANNidxArray			pidx,			// point indices
					int					n,				// number of points
					int					d,				// dimension along which to split
					ANNcoord			&cv,			// cutting value
					int					n_lo);			// split into n_lo and n-n_lo

void annPlaneSplit(				// split points by a plane
				   ANNpointArray		pa,				// points to split
				   ANNidxArray			pidx,			// point indices
				   int					n,				// number of points
				   int					d,				// dimension along which to split
				   ANNcoord			cv,				// cutting value
				   int					&br1,			// first break (values < cv)
				   int					&br2);			// second break (values == cv)

void annBoxSplit(				// split points by a box
				 ANNpointArray		pa,				// points to split
				 ANNidxArray			pidx,			// point indices
				 int					n,				// number of points
				 int					dim,			// dimension of space
				 ANNorthRect			&box,			// the box
				 int					&n_in);			// number of points inside (returned)

int annSplitBalance(			// determine balance factor of a split
					ANNpointArray		pa,				// points to split
					ANNidxArray			pidx,			// point indices
					int					n,				// number of points
					int					d,				// dimension along which to split
					ANNcoord			cv);			// cutting value

void annBox2Bnds(				// convert inner box to bounds
				 const ANNorthRect	&inner_box,		// inner box
				 const ANNorthRect	&bnd_box,		// enclosing box
				 int					dim,			// dimension of space
				 int					&n_bnds,		// number of bounds (returned)
				 ANNorthHSArray		&bnds);			// bounds array (returned)

void annBnds2Box(				// convert bounds to inner box
				 const ANNorthRect	&bnd_box,		// enclosing box
				 int					dim,			// dimension of space
				 int					n_bnds,			// number of bounds
				 ANNorthHSArray		bnds,			// bounds array
				 ANNorthRect			&inner_box);	// inner box (returned)

#endif
