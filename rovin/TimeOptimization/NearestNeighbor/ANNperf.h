#ifndef ANNperf_H
#define ANNperf_H

//----------------------------------------------------------------------
//	basic includes
//----------------------------------------------------------------------

#include "ANN.h"					// basic ANN includes


//----------------------------------------------------------------------
//		Operation count updates
//----------------------------------------------------------------------

#ifdef ANN_PERF
#define ANN_FLOP(n)	{ann_Nfloat_ops += (n);}
#define ANN_LEAF(n)	{ann_Nvisit_lfs += (n);}
#define ANN_SPL(n)	{ann_Nvisit_spl += (n);}
#define ANN_SHR(n)	{ann_Nvisit_shr += (n);}
#define ANN_PTS(n)	{ann_Nvisit_pts += (n);}
#define ANN_COORD(n)	{ann_Ncoord_hts += (n);}
#else
#define ANN_FLOP(n)
#define ANN_LEAF(n)
#define ANN_SPL(n)
#define ANN_SHR(n)
#define ANN_PTS(n)
#define ANN_COORD(n)
#endif


#endif
