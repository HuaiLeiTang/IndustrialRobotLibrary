#include "ANNx.h"					// all ANN includes
#include "ANNperf.h"				// ANN performance 

ANNdist annDist(						// interpoint squared distance
				int					dim,
				ANNpoint			p,
				ANNpoint			q)
{
	register int d;
	register ANNcoord diff;
	register ANNcoord dist;

	dist = 0;
	for (d = 0; d < dim; d++) {
		diff = p[d] - q[d];
		dist = ANN_SUM(dist, ANN_POW(diff));
	}
	ANN_FLOP(3*dim)					// performance counts
		ANN_PTS(1)
		ANN_COORD(dim)
		return dist;
}


ANNpoint annAllocPt(int dim, ANNcoord c)		// allocate 1 point
{
	ANNpoint p = new ANNcoord[dim];
	for (int i = 0; i < dim; i++) p[i] = c;
	return p;
}

ANNpointArray annAllocPts(int n, int dim)		// allocate n pts in dim
{
	ANNpointArray pa = new ANNpoint[n];			// allocate points
	ANNpoint	  p  = new ANNcoord[n*dim];		// allocate space for coords
	for (int i = 0; i < n; i++) {
		pa[i] = &(p[i*dim]);
	}
	return pa;
}

void annDeallocPt(ANNpoint &p)					// deallocate 1 point
{
	delete [] p;
	p = NULL;
}

void annDeallocPts(ANNpointArray &pa)			// deallocate points
{
	delete [] pa[0];							// dealloc coordinate storage
	delete [] pa;								// dealloc points
	pa = NULL;
}

ANNpoint annCopyPt(int dim, ANNpoint source)	// copy point
{
	ANNpoint p = new ANNcoord[dim];
	for (int i = 0; i < dim; i++) p[i] = source[i];
	return p;
}

// assign one rect to another
void annAssignRect(int dim, ANNorthRect &dest, const ANNorthRect &source)
{
	for (int i = 0; i < dim; i++) {
		dest.lo[i] = source.lo[i];
		dest.hi[i] = source.hi[i];
	}
}

// is point inside rectangle?
ANNbool ANNorthRect::inside(int dim, ANNpoint p)
{
	for (int i = 0; i < dim; i++) {
		if (p[i] < lo[i] || p[i] > hi[i]) return ANNfalse;
	}
	return ANNtrue;
}



int	ANNmaxPtsVisited = 0;	// maximum number of pts visited
int	ANNptsVisited;			// number of pts visited in search

void annMaxPtsVisit(			// set limit on max. pts to visit in search
					int					maxPts)			// the limit
{
	ANNmaxPtsVisited = maxPts;
}
