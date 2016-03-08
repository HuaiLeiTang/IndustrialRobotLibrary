#ifndef ANN_KD_SPLIT_H
#define ANN_KD_SPLIT_H

#include "kd_tree.h"					// kd-tree definitions

//----------------------------------------------------------------------
//	External entry points
//		These are all splitting procedures for kd-trees.
//----------------------------------------------------------------------

void kd_split(							// standard (optimized) kd-splitter
			  ANNpointArray		pa,				// point array (unaltered)
			  ANNidxArray			pidx,			// point indices (permuted on return)
			  const ANNorthRect	&bnds,			// bounding rectangle for cell
			  int					n,				// number of points
			  int					dim,			// dimension of space
			  int					&cut_dim,		// cutting dimension (returned)
			  ANNcoord			&cut_val,		// cutting value (returned)
			  int					&n_lo);			// num of points on low side (returned)

void midpt_split(						// midpoint kd-splitter
				 ANNpointArray		pa,				// point array (unaltered)
				 ANNidxArray			pidx,			// point indices (permuted on return)
				 const ANNorthRect	&bnds,			// bounding rectangle for cell
				 int					n,				// number of points
				 int					dim,			// dimension of space
				 int					&cut_dim,		// cutting dimension (returned)
				 ANNcoord			&cut_val,		// cutting value (returned)
				 int					&n_lo);			// num of points on low side (returned)

void sl_midpt_split(					// sliding midpoint kd-splitter
					ANNpointArray		pa,				// point array (unaltered)
					ANNidxArray			pidx,			// point indices (permuted on return)
					const ANNorthRect	&bnds,			// bounding rectangle for cell
					int					n,				// number of points
					int					dim,			// dimension of space
					int					&cut_dim,		// cutting dimension (returned)
					ANNcoord			&cut_val,		// cutting value (returned)
					int					&n_lo);			// num of points on low side (returned)

void fair_split(						// fair-split kd-splitter
				ANNpointArray		pa,				// point array (unaltered)
				ANNidxArray			pidx,			// point indices (permuted on return)
				const ANNorthRect	&bnds,			// bounding rectangle for cell
				int					n,				// number of points
				int					dim,			// dimension of space
				int					&cut_dim,		// cutting dimension (returned)
				ANNcoord			&cut_val,		// cutting value (returned)
				int					&n_lo);			// num of points on low side (returned)

void sl_fair_split(						// sliding fair-split kd-splitter
				   ANNpointArray		pa,				// point array (unaltered)
				   ANNidxArray			pidx,			// point indices (permuted on return)
				   const ANNorthRect	&bnds,			// bounding rectangle for cell
				   int					n,				// number of points
				   int					dim,			// dimension of space
				   int					&cut_dim,		// cutting dimension (returned)
				   ANNcoord			&cut_val,		// cutting value (returned)
				   int					&n_lo);			// num of points on low side (returned)

#endif
