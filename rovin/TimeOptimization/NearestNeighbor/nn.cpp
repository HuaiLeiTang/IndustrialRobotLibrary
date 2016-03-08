#include "nn.h"

int anncount = 0;

ANN::ANN(ANNpointArray points, int n1, int n2, int size, int dim, int *topology, ANNpoint scaling, int NumNeighbors) {
	int n = 0, i = 0, j = 0;

	anncount++;

	numPoints = size;
	dimension = dim;

	query_pt = annAllocPt(dimension);		// allocate query point
	data_pts = annAllocPts(numPoints, dimension); // allocate data points
	nn_idx = new ANNidx[NumNeighbors];   		// allocate near neighbor indices
	dists = new ANNdist[NumNeighbors];   		// allocate near neighbor dists
	epsilon = 0.0;

	node_indices = new int[numPoints];   		// allocate indices of the points

	// Copy nodes to data_pts
	i = 0;
	for (n = n1; n!= n2; n++) {
		node_indices[i] = n;
		for (j = 0; j < dimension; j++) {
			data_pts[i][j] = scaling[j]*points[n][j];
		}
		i++;
	}

	// Initialize ANN
	the_tree = new ANNkd_tree(data_pts,		// the data points
		numPoints,		// number of points
		dimension, 		// dimension of space
		scaling, 		// scaling of the coordinates
		topology);		// topology of the space

}

ANN::~ANN() {
	annDeallocPts(data_pts);			// deallocate data points
	annDeallocPt(query_pt);			// deallocate query points
	delete [] node_indices;			// deallocate point indices
	delete(nn_idx);				// deallocate the near neighbors indices
	delete(dists);				// deallocate the near neighbors distances
	delete(the_tree);				// deastroy the ANN tree
}


int ANN::NearestNeighbor(int *topology, double *scaling, const ANNpoint &x, double &d_best) 
{
	int i;
	int best_node;				// index of the closest node to x


	// Transfer x to query_pt
	for (i = 0; i < dimension ; i++) {
		query_pt[i] = scaling[i]*x[i];
	}

	the_tree->annkSearch(			
		query_pt,			// query point
		1,   	       			// number of near neighbors
		nn_idx,		       		// nearest neighbors (returned)
		dists,				// distance (returned)
		epsilon);	       		// error bound

	best_node = node_indices[nn_idx[0]];
	d_best = sqrt(dists[0]);
	return best_node;
}

void ANN::NearestNeighbor(int *topology, ANNpoint scaling, const ANNpoint &x, int NumNeighbors, 
						  ANNpoint &d_best, int *&n_best) 
{
	int i;  

	// Transfer x to query_pt
	for (i = 0; i < dimension ; i++) {
		query_pt[i] = scaling[i]*x[i];
	}

	the_tree->annkSearch(			
		query_pt,			// query point
		NumNeighbors,  			// number of near neighbors
		nn_idx,		       		// nearest neighbors (returned)
		dists,				// distance (returned)
		epsilon);	       		// error bound

	for (i = 0; i < NumNeighbors; i++) {
		n_best[i] = node_indices[nn_idx[i]];
		d_best[i] = sqrt(dists[i]);
	}
}

int ANN::LastNode() {
	return node_indices[numPoints-1];
}

int ANN::FirstNode() {
	return node_indices[0];
}
