#ifndef NNS_H
#define NNS_H

#define BUCKETSIZE 10

typedef struct kdtree
{
	int cutdim;
	float cutval;
	float *max, *min;
	struct kdtree *left, *right;
	int *idx, num;
} kdtree;

kdtree *kdtree_build(float **p, int num, int dim);
void kdtree_free(kdtree *t);
int kdtree_search(float *q, float **p, int dim, kdtree *t, float eps, int num, int *nlist);

#endif
