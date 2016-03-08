#ifndef ANN_kd_search_H
#define ANN_kd_search_H

#include "kd_tree.h"					// kd-tree declarations
#include "kd_util.h"					// kd-tree utilities
#include "pr_queue_k.h"					// k-element priority queue

#include "ANNperf.h"				// performance evaluation

//----------------------------------------------------------------------
//	More global variables
//		These are active for the life of each call to annkSearch(). They
//		are set to save the number of variables that need to be passed
//		among the various search procedures.
//----------------------------------------------------------------------

extern int				ANNkdDim;		// dimension of space (static copy)
extern ANNpoint			ANNkdQ;			// query point (static copy)
extern double			ANNkdMaxErr;	// max tolerable squared error
extern ANNpointArray	ANNkdPts;		// the points (static copy)
extern ANNmin_k			*ANNkdPointMK;	// set of k closest points
extern int				ANNptsVisited;	// number of points visited
extern double			*ANNScale;		// scaling array 
extern int				*ANNTopology;	// topology array 

#endif
