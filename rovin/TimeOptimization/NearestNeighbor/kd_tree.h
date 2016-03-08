#ifndef ANN_kd_tree_H
#define ANN_kd_tree_H

#include "ANNx.h"					// all ANN includes


class ANNkd_node{						// generic kd-tree node (empty shell)
public:
	virtual ~ANNkd_node() {}					// virtual distroyer

	virtual void ann_search(ANNdist) = 0;		// tree search

	friend class ANNkd_tree;					// allow kd-tree to access us
};

typedef void (*ANNkd_splitter)(			// splitting routine for kd-trees
							   ANNpointArray		pa,				// point array (unaltered)
							   ANNidxArray			pidx,			// point indices (permuted on return)
							   const ANNorthRect	&bnds,			// bounding rectangle for cell
							   int					n,				// number of points
							   int					dim,			// dimension of space
							   int					&cut_dim,		// cutting dimension (returned)
							   ANNcoord			&cut_val,		// cutting value (returned)
							   int					&n_lo);			// num of points on low side (returned)

class ANNkd_leaf: public ANNkd_node		// leaf node for kd-tree
{
	int					n_pts;			// no. points in bucket
	ANNidxArray			bkt;			// bucket of points
public:
	ANNkd_leaf(							// constructor
		int				n,				// number of points
		ANNidxArray		b)				// bucket
	{
		n_pts		= n;			// number of points in bucket
		bkt			= b;			// the bucket
	}

	~ANNkd_leaf() { }					// destructor (none)

	virtual void ann_search(ANNdist);			// standard search

};

extern ANNkd_leaf *KD_TRIVIAL;					// trivial (empty) leaf node

class ANNkd_split : public ANNkd_node	// splitting node of a kd-tree
{
	int					cut_dim;		// dim orthogonal to cutting plane
	ANNcoord			cut_val;		// location of cutting plane
	ANNcoord			cd_bnds[2];		// lower and upper bounds of
	// rectangle along cut_dim
	ANNkd_ptr			child[2];		// left and right children

	// for split in P3
	int			dim1, dim2, dim3;    	// dimensions of the rest of P3
	ANNcoord		bnds1[2];      		// bounding values of the box 
	ANNcoord		bnds2[2];      		// in these dimensions
	ANNcoord		bnds3[2];

public:
	ANNkd_split(						// constructor
		int cd,							// cutting dimension
		ANNcoord cv,					// cutting value
		ANNcoord lv, ANNcoord hv,				// low and high values
		ANNkd_ptr lc=NULL, ANNkd_ptr hc=NULL)	// children
	{
		cut_dim		= cd;					// cutting dimension
		cut_val		= cv;					// cutting value
		cd_bnds[LO] = lv;				// lower bound for rectangle
		cd_bnds[HI] = hv;				// upper bound for rectangle
		child[LO]	= lc;				// left child
		child[HI]	= hc;				// right child
	}

	ANNkd_split(				// constructor for the splitting node in P3
		int cd,					// cutting dimension
		ANNcoord cv,				// cutting value
		ANNcoord lv, ANNcoord hv,		// low and high values
		int d1, int d2, int d3,
		ANNcoord lb1, ANNcoord hb1,
		ANNcoord lb2, ANNcoord hb2,
		ANNcoord lb3, ANNcoord hb3,
		ANNkd_ptr lc=NULL, ANNkd_ptr hc=NULL)	// children
	{
		cut_dim	= cd;			// cutting dimension
		cut_val	= cv;			// cutting value
		cd_bnds[LO] = lv;			// lower bound for rectangle
		cd_bnds[HI] = hv;			// upper bound for rectangle
		child[LO]	= lc;			// left child
		child[HI]	= hc;			// right child
		dim1 = d1;				// rest of the dimensions of P3
		dim2 = d2;
		dim3 = d3;
		bnds1[LO] = lb1;			// bounding values in d1 and d2
		bnds1[HI] = hb1;
		bnds2[LO] = lb2;
		bnds2[HI] = hb2;
		bnds3[LO] = lb3;
		bnds3[HI] = hb3;
	}

	~ANNkd_split()				// destructor
	{
		if (child[LO]!= NULL && child[LO]!= KD_TRIVIAL) delete child[LO];
		if (child[HI]!= NULL && child[HI]!= KD_TRIVIAL) delete child[HI];
	}

	virtual void ann_search(ANNdist);			// standard search
};

ANNkd_ptr rkd_tree(				// recursive construction of kd-tree
				   ANNpointArray		pa,				// point array (unaltered)
				   ANNidxArray			pidx,			// point indices to store in subtree
				   int					n,				// number of points
				   int					dim,			// dimension of space
				   int					bsp,			// bucket space
				   ANNorthRect			&bnd_box,		// bounding box for current node
				   ANNkd_splitter		splitter);		// splitting routine

#endif
