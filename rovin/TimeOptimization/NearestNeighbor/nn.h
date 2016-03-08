#ifndef DNN_NN_H
#define DNN_NN_H

#include "ANN.h"

class ANN {

protected:
	int numPoints;  				// Total number currently stored in the ANN tree
	int dimension;				// dimension of the space
	ANNpointArray	data_pts;			// data points
	ANNpoint query_pt;		        	// query point
	ANNidxArray nn_idx;				// near neighbor indices
	ANNdistArray dists;				// near neighbor distances
	ANNkd_tree *the_tree;		        	// search structure
	int *node_indices;        			// indices to the points in MultiANN
	double epsilon;        			// the error bound for the search

public:
	ANN() {};
	ANN(						// constructor
		ANNpointArray points, 			// array of data points
		int n1, 					// starting index in points	
		int n2, 					// ending index in points
		int size, 				// size of the chunk of points to be used
		int dim, 					// dimension of the space
		int *topology, 				// topology of the space
		ANNpoint scaling, 			// scaling of the coordinates
		int NumNeighbors);			// number of nearest neighbors to be returned

	~ANN();					// destructor

	int  NearestNeighbor(  			// 1-nearest neighbor query
		int *topology, 		// topology of the space
		ANNpoint scaling, 	// scaling of the coordinates
		const ANNpoint &x,	// query point 
		double &d_best); 	// distance from the nearest neighbor to x (returned)

	void NearestNeighbor(  			// k-nearest neighbor query
		int *topology,  		// topology of the space
		ANNpoint scaling, 	// scaling of the coordinates 
		const ANNpoint &x, 	// query point 
		int NumNeighbors, 	// the number of nearest neighbors to be returned (k)
		ANNpoint &d_best, 	// array of distances from the k nearest neighbors to x (returned) 
		int *&n_best); 		// array of k nearest neighbors to x (returned) 

	int LastNode();  				// The index of the last node in MultiANN that is stored in the ANN
	int FirstNode(); 				// The index of the first node in MultiANN that is stored in the ANN
};
#endif

