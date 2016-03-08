#ifndef PR_QUEUE_K_H
#define PR_QUEUE_K_H

#include "ANNx.h"					// all ANN includes
#include "ANNperf.h"				// performance evaluation

typedef ANNdist			PQKkey;			// key field is distance
typedef int				PQKinfo;		// info field is int


const PQKkey	PQ_NULL_KEY  =  ANN_DIST_INF;	// nonexistent key value
const PQKinfo	PQ_NULL_INFO =  ANN_NULL_IDX;	// nonexistent info value

class ANNmin_k {
	struct mk_node {					// node in min_k structure
		PQKkey			key;			// key value
		PQKinfo			info;			// info field (user defined)
	};

	int			k;						// max number of keys to store
	int			n;						// number of keys currently active
	mk_node		*mk;					// the list itself

public:
	ANNmin_k(int max)					// constructor (given max size)
	{
		n = 0;						// initially no items
		k = max;					// maximum number of items
		mk = new mk_node[max+1];	// sorted array of keys
	}

	~ANNmin_k()							// destructor
	{ delete [] mk; }

	PQKkey ANNmin_key()					// return minimum key
	{ return (n > 0 ? mk[0].key : PQ_NULL_KEY); }

	PQKkey max_key()					// return maximum key
	{ return (n == k ? mk[k-1].key : PQ_NULL_KEY); }

	PQKkey ith_smallest_key(int i)		// ith smallest key (i in [0..n-1])
	{ return (i < n ? mk[i].key : PQ_NULL_KEY); }

	PQKinfo ith_smallest_info(int i)	// info for ith smallest (i in [0..n-1])
	{ return (i < n ? mk[i].info : PQ_NULL_INFO); }

	inline void insert(					// insert item (inlined for speed)
		PQKkey kv,						// key value
		PQKinfo inf)					// item info
	{
		register int i;
		// slide larger values up
		for (i = n; i > 0; i--) {
			if (mk[i-1].key > kv)
				mk[i] = mk[i-1];
			else
				break;
		}
		mk[i].key = kv;				// store element here
		mk[i].info = inf;
		if (n < k) n++;				// increment number of items
		ANN_FLOP(k-i+1)				// increment floating ops
	}
};

#endif
