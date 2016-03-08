#include "kd_search.h"					// kd-search declarations

#define PI 3.1415926535897932385


int				ANNkdDim;				// dimension of space
ANNpoint		ANNkdQ;					// query point
double			ANNkdMaxErr;			// max tolerable squared error
ANNpointArray	ANNkdPts;				// the points
ANNmin_k		*ANNkdPointMK;			// set of k closest points
double			*ANNScale;				// scaling array
int				*ANNTopology;			// topology array

//----------------------------------------------------------------------
//	annkSearch - search for the k nearest neighbors
//----------------------------------------------------------------------

void ANNkd_tree::annkSearch(
							ANNpoint			q,				// the query point
							int					k,				// number of near neighbors to return
							ANNidxArray			nn_idx,			// nearest neighbor indices (returned)
							ANNdistArray		dd,				// the approximate nearest neighbor
							double				eps)			// the error bound
{

	ANNkdDim = dim;			// copy arguments to static equivs
	ANNScale = new double [dim];
	ANNTopology = new int [dim];
	for (int i = 0; i < dim; i++){
		ANNScale[i] = Scale[i];
		ANNTopology[i] = Topology[i];
	}
	ANNkdQ = q;
	ANNkdPts = pts;
	ANNptsVisited = 0;			// initialize count of points visited

	ANNkdMaxErr = ANN_POW(1.0 + eps);
	//    FLOP(2)				// increment floating op count

	ANNkdPointMK = new ANNmin_k(k);	// create set for closest k points
	// search starting at the root
	root->ann_search(annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim, ANNScale, ANNTopology));

	for (int i = 0; i < k; i++) {	// extract the k-th closest points
		dd[i] = ANNkdPointMK->ith_smallest_key(i);
		nn_idx[i] = ANNkdPointMK->ith_smallest_info(i);
	}
	delete [] ANNScale;
	delete [] ANNTopology;
	delete ANNkdPointMK;		// deallocate closest point set

}

//----------------------------------------------------------------------
//	kd_split::ann_search - search a splitting node
//----------------------------------------------------------------------

void ANNkd_split::ann_search(ANNdist box_dist)
{

	// check dist calc termination condition
	if (ANNmaxPtsVisited && ANNptsVisited > ANNmaxPtsVisited) return;

	if (ANNTopology[cut_dim] == 1){
		// distance to cutting plane
		ANNcoord cut_diff = ANNkdQ[cut_dim] - cut_val;

		if (cut_diff < 0) {			// left of cutting plane
			child[LO]->ann_search(box_dist);// visit closer child first

			ANNcoord box_diff = cd_bnds[LO] - ANNkdQ[cut_dim];
			if (box_diff < 0)		// within bounds - ignore
				box_diff = 0;
			// distance to further box
			box_dist = (ANNdist) ANN_SUM(box_dist,
				ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

			// visit further child if close enough
			if (box_dist * ANNkdMaxErr < ANNkdPointMK->max_key())
				child[HI]->ann_search(box_dist);

		}
		else {				// right of cutting plane
			child[HI]->ann_search(box_dist);// visit closer child first

			ANNcoord box_diff = ANNkdQ[cut_dim] - cd_bnds[HI];
			if (box_diff < 0)		// within bounds - ignore
				box_diff = 0;
			// distance to further box
			box_dist = (ANNdist) ANN_SUM(box_dist,
				ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

			// visit further child if close enough
			if (box_dist * ANNkdMaxErr < ANNkdPointMK->max_key())
				child[LO]->ann_search(box_dist);

		}
		//FLOP(10)				// increment floating ops
		//SPL(1)				// one more splitting node visited
	}
	else if (ANNTopology[cut_dim] == 2) {
		ANNdist box_dist1, box_dist2;

		// distance to cutting plane
		ANNcoord cut_diff1 = ANNkdQ[cut_dim] - cut_val;

		if (cut_diff1 < 0) {       		// left of cutting plane

			ANNcoord cut_diff2 = 2*PI*ANNScale[cut_dim] + (ANNkdQ[cut_dim] - cut_val);

			ANNcoord box_diff1 = cd_bnds[LO] - ANNkdQ[cut_dim];
			if (box_diff1 < 0)		// within bounds - ignore
				box_diff1 = 0;

			ANNcoord box_diff2 = 2*PI*ANNScale[cut_dim] - (cd_bnds[HI] - ANNkdQ[cut_dim]);
			if (box_diff2 < 0)		// within bounds - ignore
				box_diff2 = 0;

			if (box_diff1 < box_diff2) {
				child[LO]->ann_search(box_dist);// visit closer child first
				// distance to further box
				box_dist1 = (ANNdist) ANN_SUM(box_dist,
					ANN_DIFF(ANN_POW(box_diff1), ANN_POW(cut_diff1)));


				// distance to further box
				box_dist2 = (ANNdist) ANN_SUM(box_dist,
					ANN_DIFF(ANN_POW(box_diff1), ANN_POW(box_diff2)));
				if (box_dist1 < box_dist2) {

					// visit further child if close enough
					if (box_dist1 * ANNkdMaxErr < ANNkdPointMK->max_key())
						child[HI]->ann_search(box_dist1);
				}
				else {
					if (box_dist2 * ANNkdMaxErr < ANNkdPointMK->max_key())
						child[HI]->ann_search(box_dist2);
				}
			}
			else {
				child[HI]->ann_search(box_dist);// visit closer child first
				// distance to further box
				box_dist1 = (ANNdist) ANN_SUM(box_dist,
					ANN_DIFF(ANN_POW(box_diff2), ANN_POW(cut_diff2)));

				// distance to further box
				box_dist2 = (ANNdist) ANN_SUM(box_dist,
					ANN_DIFF(ANN_POW(box_diff2), ANN_POW(box_diff1)));
				if (box_dist1 < box_dist2) {

					// visit further child if close enough
					if (box_dist1 * ANNkdMaxErr < ANNkdPointMK->max_key())
						child[LO]->ann_search(box_dist1);
				}
				else {
					if (box_dist2 * ANNkdMaxErr < ANNkdPointMK->max_key())
						child[LO]->ann_search(box_dist2);
				}

			}
		}
		else {				// right of cutting plane
			ANNcoord cut_diff2 = 2*PI*ANNScale[cut_dim] - (ANNkdQ[cut_dim] - cut_val);

			ANNcoord box_diff1 =  ANNkdQ[cut_dim] - cd_bnds[HI];
			if (box_diff1 < 0)		// within bounds - ignore
				box_diff1 = 0;

			ANNcoord box_diff2 = 2*PI*ANNScale[cut_dim] - (ANNkdQ[cut_dim] - cd_bnds[LO]);
			if (box_diff2 < 0)		// within bounds - ignore
				box_diff2 = 0;

			if (box_diff1 < box_diff2) {
				child[HI]->ann_search(box_dist);// visit closer child first
				// distance to further box
				box_dist1 = (ANNdist) ANN_SUM(box_dist,
					ANN_DIFF(ANN_POW(box_diff1), ANN_POW(cut_diff1)));


				// distance to further box
				box_dist2 = (ANNdist) ANN_SUM(box_dist,
					ANN_DIFF(ANN_POW(box_diff1), ANN_POW(box_diff2)));
				if (box_dist1 < box_dist2) {

					// visit further child if close enough
					if (box_dist1 * ANNkdMaxErr < ANNkdPointMK->max_key())
						child[LO]->ann_search(box_dist1);
				}
				else {
					if (box_dist2 * ANNkdMaxErr < ANNkdPointMK->max_key())
						child[LO]->ann_search(box_dist2);
				}
			}
			else {
				child[LO]->ann_search(box_dist);// visit closer child first
				// distance to further box
				box_dist1 = (ANNdist) ANN_SUM(box_dist,
					ANN_DIFF(ANN_POW(box_diff2), ANN_POW(cut_diff2)));


				// distance to further box
				box_dist2 = (ANNdist) ANN_SUM(box_dist,
					ANN_DIFF(ANN_POW(box_diff2), ANN_POW(box_diff1)));
				if (box_dist1 < box_dist2) {

					// visit further child if close enough
					if (box_dist1 * ANNkdMaxErr < ANNkdPointMK->max_key())
						child[HI]->ann_search(box_dist1);
				}
				else {
					if (box_dist2 * ANNkdMaxErr < ANNkdPointMK->max_key())
						child[HI]->ann_search(box_dist2);
				}

			}
		}
		//FLOP(10)				// increment floating ops
		//SPL(1)				// one more splitting node visited
	}

	else {

		ANNdist box_dist1, box_dist2;

		ANNRectangle *RectLO = new ANNRectangle(cd_bnds[LO], cut_val, bnds1[LO], bnds1[HI], 
			bnds2[LO], bnds2[HI], bnds3[LO], bnds3[HI]);
		ANNRectangle *RectHI = new ANNRectangle(cut_val, cd_bnds[HI], bnds1[LO], bnds1[HI], 
			bnds2[LO], bnds2[HI], bnds3[LO], bnds3[HI]);
		ANNdist distLO, distLO1, distHI, distHI1;
		distLO = RectLO->DistPointRectangle (ANNkdQ[cut_dim], ANNkdQ[dim1], ANNkdQ[dim2], ANNkdQ[dim3]);
		distLO1 = RectLO->DistPointRectangle (-ANNkdQ[cut_dim], -ANNkdQ[dim1], -ANNkdQ[dim2], -ANNkdQ[dim3]);
		if (distLO1 < distLO) distLO = distLO1;

		distHI = RectHI->DistPointRectangle (ANNkdQ[cut_dim], ANNkdQ[dim1], ANNkdQ[dim2], ANNkdQ[dim3]);
		distHI1 = RectHI->DistPointRectangle (-ANNkdQ[cut_dim], -ANNkdQ[dim1], -ANNkdQ[dim2], -ANNkdQ[dim3]);
		if (distHI1 < distHI) distHI = distHI1;

		if (distLO < distHI) {
			box_dist1 = box_dist;
			box_dist2 = box_dist - distLO + distHI;
			child[LO]->ann_search(box_dist1);

			if (box_dist2 * ANNkdMaxErr < ANNkdPointMK->max_key())
				child[HI]->ann_search(box_dist2);
		}
		else {
			box_dist1 = box_dist;
			box_dist2 = box_dist - distHI + distLO;
			child[HI]->ann_search(box_dist1);

			if (box_dist2 * ANNkdMaxErr < ANNkdPointMK->max_key())
				child[LO]->ann_search(box_dist2);
		}
		delete RectLO;
		delete RectHI;
		//FLOP(10)				// increment floating ops
		//SPL(1)				// one more splitting node visited

	}
}


//----------------------------------------------------------------------
//	kd_leaf::ann_search - search points in a leaf node
//		Note: The unreadability of this code is the result of
//		some fine tuning to replace indexing by pointer operations.
//----------------------------------------------------------------------

void ANNkd_leaf::ann_search(ANNdist box_dist)
{

	register ANNdist dist;		// distance to data point
	register ANNcoord* pp;		// data coordinate pointer
	register ANNcoord* qq;		// query coordinate pointer
	register ANNdist min_dist;		// distance to k-th closest point
	register int d;

	register ANNcoord t, t1;
	double cosTheta, dotProduct, norm1, norm2;

	min_dist = ANNkdPointMK->max_key();	// k-th smallest distance so far

	for (int i = 0; i < n_pts; i++) {	// check points in bucket

		pp = ANNkdPts[bkt[i]];		// first coord of next data point
		qq = ANNkdQ; 	    		// first coord of query point
		dist = 0;

		for(d = 0; d < ANNkdDim; d++) {
			//COORD(1)			// one more coordinate hit
			//FLOP(4)			// increment floating ops

			if (ANNTopology[d] == 2) {
				t = fabs(*(qq) - *(pp));    	// compute length and adv coordinate
				t1 = 2*PI*ANNScale[d] - t;
				if (t > t1) t = t1;
				t = ANN_POW(t);
			}
			else if (ANNTopology[d] == 1){
				t = ANN_POW(*(qq) - *(pp));	// compute length and adv coordinate
			}
			else {
				dotProduct = *(qq) * *(pp) + *(qq+1) * *(pp+1) + *(qq+2) * *(pp+2) + *(qq+3) * *(pp+3);
				cosTheta = dotProduct/ANN_POW(ANNScale[d]);

				if (cosTheta > 1) {
					norm1 = *(qq) * *(qq) + *(qq+1) * *(qq+1) + *(qq+2) * *(qq+2) + *(qq+3) * *(qq+3);
					norm2 = *(pp) * *(pp) + *(pp+1) * *(pp+1) + *(pp+2) * *(pp+2) + *(pp+3) * *(pp+3);
					cosTheta = dotProduct/(sqrt(norm1*norm2));
				}
				t = acos(cosTheta);
				t1 = acos(-cosTheta);
				if (t > t1) t = t1;
				t = ANN_POW(ANNScale[d]*t);
				d = d + 3; qq++;qq++;qq++; pp++;pp++;pp++;
			}

			dist = ANN_SUM(dist, t);
			qq++; pp++;

			// exceeds dist to k-th smallest?
			if( dist > min_dist) {
				break;
			}
		}

		if (d >= ANNkdDim &&			// among the k best?
			(ANN_ALLOW_SELF_MATCH || dist!=0)) {	// and no self-match problem
				// add it to the list
				ANNkdPointMK->insert(dist, bkt[i]);
				min_dist = ANNkdPointMK->max_key();
		}
	}
	//LEAF(1)				// one more leaf node visited
	//PTS(n_pts)				// increment points visited
	ANNptsVisited += n_pts;		// increment number of points visited

}
