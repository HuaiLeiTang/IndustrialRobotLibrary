#ifndef DNN_MULTIANN_H
#define DNN_MULTIANN_H

#include <cmath>

#include "nn.h"

#define ANN_STARTING_INDEX 4  			// Start a new ANN for each 2^n samples
#define ANN_MAXIMUM_INDEX 30  			// Absolute max on number of nodes (2^n)

class MultiANN{
public:

	int size;					// the number of points in MultiANN
	int dimension;				// dimension of the space
	int LastNodeCount;				// the index of the last node
	ANNpointArray points_coor;			// array of the coordinates of data points
	int *points_ptr;				// array of the pointers to data points
	int NumNeighbors;				// number of nearest neighbors to be returned (k) 
	int *topology;				// topology of the space 
	ANNpoint scaling;				// scaling of the coordinates

	ANN  *AnnArray[ANN_MAXIMUM_INDEX+1];		// array of ANNs to hold data points

	MultiANN() {};				// constructor
	MultiANN(int dim, int Maxpoint);										// Âù¼ö ¸¸µë
	MultiANN(int dim, int Maxpoint, ANNpoint x_coor, int *x_ptr);         // Âù¼ö ¸¸µë
	MultiANN(					// constructor
		int dim, 				// dimension of the space
		ANNpoint x_coor, 		       	// coordinate of the initial point in the data structure
		int *x_ptr, 		       	// pointer to the initial point in the data structure
		int k, 				// number of nearest neighbors to be returned
		int *topology, 			// topology of the space
		ANNpoint scaling);			// scaling of the coordinates

	MultiANN(					// constructor
		int dim,  				// dimension of the space
		int k, 				// number of nearest neighbors to be returned 
		int *topology, 			// topology of the space 
		ANNpoint scaling);			// scaling of the coordinates

	~MultiANN();					// destructor
	void ResetMultiANN();				// destroys all the arrays of data points
	void AddPoint(				// dynamic update of the data points
		ANNpoint x_coor,		// the point's coordinates 
		int *x_ptr);			// the pointer to the point 

	void UpdateAnnArray();			// updates the arrays of data points

	void NearestNeighbor(
	//int NearestNeighbor(  			// 1-nearest neighbor query
		const ANNpoint &x,	// query point  
		int    &best_idx, 	// distance from the nearest neighbor to x (returned)
		double* best_dist); // distance from the nearest neighbor to x (returned)

	void NearestNeighbor(  			// k-nearest neighbor query
		const ANNpoint &x, 	// query point 
		ANNpoint &best_dist,  	// array of distances from the k nearest neighbors to x (returned) 
		int *&best_idx, 		// array of indices of k nearest neighbors to x (returned) 
		int *&best_ptr);      	// array of pointers to k nearest neighbors to x (returned) 
};
#endif
