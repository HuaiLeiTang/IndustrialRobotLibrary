#include "kd_tree.h"					// kd-tree declarations
#include "kd_split.h"					// kd-tree splitting rules
#include "kd_util.h"					// kd-tree utilities
#include "ANNperf.h"				// performance evaluation

static int				IDX_TRIVIAL[] = {0};	// trivial point index
ANNkd_leaf				*KD_TRIVIAL = NULL;		// trivial leaf node
int 					*TreeTopology;
int 					*TreeP3Topology;

#define MAX(a,b)		((a) > (b) ? (a) : (b))

const double ANN_AR_TOOBIG = 1000;				// too big an aspect ratio

ANNkd_tree::~ANNkd_tree()				// tree destructor
{
	if (root != NULL) delete root;
	if (pidx != NULL) delete [] pidx;
	if (bnd_box_lo != NULL) annDeallocPt(bnd_box_lo);
	if (bnd_box_hi != NULL) annDeallocPt(bnd_box_hi);
}

void annClose()				// close use of ANN
{
	if (KD_TRIVIAL != NULL) {
		delete KD_TRIVIAL;
		KD_TRIVIAL = NULL;
	}
}

void ANNkd_tree::SkeletonTree(			// construct skeleton tree
							  int n,							// number of points
							  int dd,							// dimension
							  int bs,							// bucket size
							  ANNpointArray pa,				// point array
							  ANNidxArray pi)					// point indices
{
	dim = dd;							// initialize basic elements
	n_pts = n;
	bkt_size = bs;
	pts = pa;							// initialize points array

	root = NULL;						// no associated tree yet

	if (pi == NULL) {					// point indices provided?
		pidx = new ANNidx[n];			// no, allocate space for point indices
		for (int i = 0; i < n; i++) {
			pidx[i] = i;				// initially identity
		}
	}
	else {
		pidx = pi;						// yes, use them
	}

	bnd_box_lo = bnd_box_hi = NULL;		// bounding box is nonexistent
	if (KD_TRIVIAL == NULL)				// no trivial leaf node yet?
		KD_TRIVIAL = new ANNkd_leaf(0, IDX_TRIVIAL);	// allocate it
}

ANNkd_tree::ANNkd_tree(					// basic constructor
					   int n,							// number of points
					   int dd,							// dimension
					   int bs)							// bucket size
{  SkeletonTree(n, dd, bs);  }			// construct skeleton tree

ANNkd_ptr rkd_tree(				// recursive construction of kd-tree
				   ANNpointArray		pa,				// point array
				   ANNidxArray			pidx,			// point indices to store in subtree
				   int					n,				// number of points
				   int					dim,			// dimension of space
				   int					bsp,			// bucket space
				   ANNorthRect			&bnd_box,		// bounding box for current node
				   ANNkd_splitter		splitter)		// splitting routine
{
	if (n <= bsp) {						// n small, make a leaf node
		if (n == 0)						// empty leaf node
			return KD_TRIVIAL;			// return (canonical) empty leaf
		else							// construct the node and return
			return new ANNkd_leaf(n, pidx); 
	}
	else {								// n large, make a splitting node
		int cd;							// cutting dimension
		ANNcoord cv;					// cutting value
		int n_lo;						// number on low side of cut
		ANNkd_node *lo, *hi;			// low and high children

		// invoke splitting procedure
		(*splitter)(pa, pidx, bnd_box, n, dim, cd, cv, n_lo);

		ANNcoord lv = bnd_box.lo[cd];	// save bounds for cutting dimension
		ANNcoord hv = bnd_box.hi[cd];

		bnd_box.hi[cd] = cv;			// modify bounds for left subtree
		lo = rkd_tree(					// build left subtree
			pa, pidx, n_lo,			// ...from pidx[0..n_lo-1]
			dim, bsp, bnd_box, splitter);
		bnd_box.hi[cd] = hv;			// restore bounds

		bnd_box.lo[cd] = cv;			// modify bounds for right subtree
		hi = rkd_tree(					// build right subtree
			pa, pidx + n_lo, n-n_lo,// ...from pidx[n_lo..n-1]
			dim, bsp, bnd_box, splitter);
		bnd_box.lo[cd] = lv;			// restore bounds

		// create the splitting node
		//ANNkd_split *ptr;
		ANNkd_split *ptr = NULL;
		if (TreeTopology[cd] == 3) {
			if (TreeP3Topology[cd] == 0)
				ptr = new ANNkd_split(cd, cv, lv, hv, cd + 1, cd + 2, cd + 3,
				bnd_box.lo[cd + 1], bnd_box.hi[cd + 1],
				bnd_box.lo[cd + 2], bnd_box.hi[cd + 2],
				bnd_box.lo[cd + 3], bnd_box.hi[cd + 3],
				lo, hi);  
			if (TreeP3Topology[cd] == 1)
				ptr = new ANNkd_split(cd, cv, lv, hv, cd - 1, cd + 1, cd + 2,
				bnd_box.lo[cd - 1], bnd_box.hi[cd - 1],
				bnd_box.lo[cd + 1], bnd_box.hi[cd + 1],
				bnd_box.lo[cd + 2], bnd_box.hi[cd + 2],
				lo, hi);  
			if (TreeP3Topology[cd] == 2)
				ptr = new ANNkd_split(cd, cv, lv, hv, cd - 2, cd - 1, cd + 1,
				bnd_box.lo[cd - 2], bnd_box.hi[cd - 2],
				bnd_box.lo[cd - 1], bnd_box.hi[cd - 1],
				bnd_box.lo[cd + 1], bnd_box.hi[cd + 1],
				lo, hi);  
			if (TreeP3Topology[cd] == 3)
				ptr = new ANNkd_split(cd, cv, lv, hv, cd - 3, cd - 2, cd - 1,
				bnd_box.lo[cd - 3], bnd_box.hi[cd - 3],
				bnd_box.lo[cd - 2], bnd_box.hi[cd - 2],
				bnd_box.lo[cd - 1], bnd_box.hi[cd - 1],
				lo, hi);  
		}
		else {
			ptr = new ANNkd_split(cd, cv, lv, hv, lo, hi);
		}		

		return ptr;						// return pointer to this node
	}
} 

ANNkd_tree::ANNkd_tree(					// construct from point array
					   ANNpointArray		pa,				// point array (with at least n pts)
					   int					n,				// number of points
					   int					dd,				// dimension
					   double				*ss,			// scaling 
					   int			 		*tt,			// topology of space
					   int					bs,				// bucket size
					   ANNsplitRule		split)			// splitting method
{
	SkeletonTree(n, dd, bs);			// set up the basic stuff
	pts = pa;							// where the points are
	if (n == 0) return;					// no points--no sweat

	Scale = new double [dim];
	Topology = new int [dim];
	TreeTopology = new int [dim];
	TreeP3Topology = new int [dim];
	for (int i = 0; i < dim; i++) {
		Scale[i] = ss[i];
		Topology[i] = tt[i];
		TreeTopology[i] = tt[i];
	}
	for (int i = 0; i < dim; i++) {
		if (tt[i] != 3) TreeP3Topology[i] = 0;
		else {
			TreeP3Topology[i++] = 0;
			TreeP3Topology[i++] = 1;
			TreeP3Topology[i++] = 2;
			TreeP3Topology[i] = 3;
		}
	}

	ANNorthRect bnd_box(dd);			// bounding box for points
	annEnclRect(pa, pidx, n, dd, bnd_box);// construct bounding rectangle
	// copy to tree structure
	bnd_box_lo = annCopyPt(dd, bnd_box.lo);
	bnd_box_hi = annCopyPt(dd, bnd_box.hi);

	switch (split) {					// build by rule
	case ANN_KD_STD:					// standard kd-splitting rule
		root = rkd_tree(pa, pidx, n, dd, bs, bnd_box, kd_split);
		break;
	case ANN_KD_MIDPT:					// midpoint split
		root = rkd_tree(pa, pidx, n, dd, bs, bnd_box, midpt_split);
		break;
	case ANN_KD_FAIR:					// fair split
		root = rkd_tree(pa, pidx, n, dd, bs, bnd_box, fair_split);
		break;
	case ANN_KD_SUGGEST:				// best (in our opinion)
	case ANN_KD_SL_MIDPT:				// sliding midpoint split
		root = rkd_tree(pa, pidx, n, dd, bs, bnd_box, sl_midpt_split);
		break;
	case ANN_KD_SL_FAIR:				// sliding fair split
		root = rkd_tree(pa, pidx, n, dd, bs, bnd_box, sl_fair_split);
		break;

	}

	delete [] TreeTopology;
	delete [] TreeP3Topology;

}
