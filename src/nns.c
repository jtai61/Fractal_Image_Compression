#include "base.h"

static heapT *h = NULL;
static int hsize;
static int cutdim;
static float **p;

/* Heap-functions for kdtree_search */

static void hinit()
{
	if (h == NULL)
		h = (heapT *)malloc(sizeof(heapT) * HEAPMAX);
	hsize = 0;
}

static void hfree()
{
	if (h != NULL)
		free(h);
}

static int hempty()
{
	return (hsize == 0);
}

static void hpush(kdtree *t, float *q, int dim)
{
	int i, pos;
	float d, tmp;

	assert(hsize < HEAPMAX);
	d = 0.0;
	for (i = 0; i < dim; i++)
	{
		if (q[i] < t->min[i])
		{
			tmp = t->min[i] - q[i];
			d += tmp * tmp;
		}
		else if (q[i] > t->max[i])
		{
			tmp = q[i] - t->max[i];
			d += tmp * tmp;
		}
	}
	pos = ++hsize;
	while (pos != 1 && h[pos / 2].dist > d)
	{
		h[pos] = h[pos / 2];
		pos /= 2;
	}
	h[pos].dist = d;
	h[pos].t = t;
}

static kdtree *hpop()
{
	int pos;
	kdtree *tmp;

	assert(hsize != 0);
	tmp = h[1].t;
	if (--hsize != 0)
	{
		pos = 1;
		while (2 * pos <= hsize)
		{
			if (2 * pos + 1 <= hsize && h[2 * pos + 1].dist < h[2 * pos].dist)
			{
				if (h[2 * pos + 1].dist < h[hsize + 1].dist)
				{
					h[pos] = h[2 * pos + 1];
					pos = 2 * pos + 1;
					continue;
				}
				else
					break;
			}
			if (h[2 * pos].dist < h[hsize + 1].dist)
			{
				h[pos] = h[2 * pos];
				pos *= 2;
				continue;
			}
			break;
		}
		h[pos] = h[hsize + 1];
	}
	return tmp;
}

/* Find 'num' (1+'eps')-nearest neighbors in kdtree 't' (built from points 'p', and put them in 'nlist') */
int kdtree_search(float *q, float **p, int dim, kdtree *t, float eps, int num, int *nlist)
{
	int i, j, found;
	float d, tmp, *ptmp, *qtmp, eps2;
	float dist2[MAX_NEIGHBOURS];

	hinit();
	hpush(t, q, dim);
	found = 0;
	eps2 = (1.0 + eps) * (1.0 + eps);
	while (!hempty() && (found < num || eps2 * h[1].dist < dist2[num - 1]))
	{
		t = hpop();
		while (t->num == 0)
		{
			if (q[t->cutdim] < t->cutval)
			{
				hpush(t->right, q, dim);
				t = t->left;
			}
			else
			{
				hpush(t->left, q, dim);
				t = t->right;
			}
		}

		/* Leaf node: process members */
		for (i = 0; i < t->num; i++)
		{
			d = 0.0;
			for (qtmp = q, ptmp = p[t->idx[i]], j = 0; j < dim; j++)
			{
				tmp = (*qtmp++) - (*ptmp++);
				d += tmp * tmp;
			}
			if (found == num && d > dist2[num - 1])
				continue;
			if (found < num)
				found++;
			for (j = found - 1; j > 0 && dist2[j - 1] > d; j--)
			{
				dist2[j] = dist2[j - 1];
				nlist[j] = nlist[j - 1];
			}
			dist2[j] = d;
			nlist[j] = t->idx[i];
		}
	}
	return found;
}

/* Functions for building a kd-tree */

static int compare(const void *a, const void *b)
{
	float p1, p2;

	p1 = p[*((int *)a)][cutdim];
	p2 = p[*((int *)b)][cutdim];
	return p1 < p2 ? -1 : p1 > p2;
}

static kdtree *kdtree_build2(int *idx, int num, int dim)
{
	int i, d;
	float spread, min, max;
	kdtree *tmp;

	if (num == 0)
		return NULL;

	tmp = (kdtree *)malloc(sizeof(kdtree));

	/* 1. Find dimension with greatest spread */

	tmp->min = (float *)malloc(sizeof(float) * dim);
	tmp->max = (float *)malloc(sizeof(float) * dim);

	spread = -1.0;
	for (d = 0; d < dim; d++)
	{
		min = max = p[idx[0]][d];
		for (i = 1; i < num; i++)
		{
			if (min > p[idx[i]][d])
				min = p[idx[i]][d];
			else if (max < p[idx[i]][d])
				max = p[idx[i]][d];
		}
		if (max - min > spread)
		{
			spread = max - min;
			cutdim = d;
		}
		tmp->min[d] = min;
		tmp->max[d] = max;
	}

	/* 2a. If the number of points is less then BUCKETSIZE:
		  pack all points in the bucket, and we're done */

	if (num <= BUCKETSIZE)
	{
		tmp->cutdim = -1;
		tmp->left = tmp->right = NULL;
		tmp->idx = (int *)malloc(sizeof(int) * num);
		tmp->num = num;
		for (i = 0; i < num; i++)
			tmp->idx[i] = idx[i];
		return tmp;
	}

	/* 2b. Otherwise, we have to build a tree recursively */

	/* 3. Sort vertices according to coordinate 'cutdim' */

	qsort(idx, num, sizeof(int), compare);

	/* 4. Split vertices evenly and recursively */

	tmp->cutdim = cutdim;
	tmp->cutval = (p[idx[num / 2 - 1]][cutdim] + p[idx[num / 2]][cutdim]) / 2.0;
	tmp->idx = NULL;
	tmp->num = 0;

	tmp->left = kdtree_build2(idx, num / 2, dim);
	tmp->right = kdtree_build2(idx + num / 2, num - num / 2, dim);

	return tmp;
}

kdtree *kdtree_build(float **p2, int num, int dim)
{
	int *idx;
	int i;
	kdtree *tmp;

	p = p2;
	idx = (int *)malloc(sizeof(int) * num);
	for (i = 0; i < num; i++)
		idx[i] = i;
	tmp = kdtree_build2(idx, num, dim);
	free(idx);
	return tmp;
}

/* To free the memory allocated by a kd-tree */

void kdtree_free(kdtree *t)
{
	if (t == NULL)
		return;
	kdtree_free(t->left);
	kdtree_free(t->right);
	if (t->idx)
		free(t->idx);
	free(t->max);
	free(t->min);
	free(t);
	hfree();
}

double manhattan_dist(double *x1, double *x2, int size)
{
	double d = 0;

	for (int i = 0; i < size; i++)
		d += fabs(x1[i] - x2[i]);
	return (d);
}

double min_dist(t_btree *tree, int i_node, double *pt)
{
	double dist_pt;

	dist_pt = manhattan_dist(pt, tree->node_bounds[0][i_node], tree->n_features);
	return (fmax(0.0, dist_pt - tree->node_data[i_node].radius));
}

void dual_swap(double *darr, int *iarr, int i1, int i2)
{
	double dtmp;
	int itmp;

	dtmp = darr[i1];
	darr[i1] = darr[i2];
	darr[i2] = dtmp;
	itmp = iarr[i1];
	iarr[i1] = iarr[i2];
	iarr[i2] = itmp;
}

t_nheap *nheap_init(int n_pts, int n_nbrs)
{
	t_nheap *h;
	int i, j;

	h = (t_nheap *)malloc(sizeof(t_nheap));
	h->n_pts = n_pts;
	h->n_nbrs = n_nbrs;
	h->distances = (double **)malloc(sizeof(double *) * n_pts);
	for (i = 0; i < n_pts; i++)
	{
		h->distances[i] = (double *)malloc(sizeof(double) * n_nbrs);
		for (j = 0; j < n_nbrs; j++)
			h->distances[i][j] = INFINITY;
	}
	h->indices = (int **)malloc(sizeof(int *) * n_pts);
	for (i = 0; i < n_pts; i++)
		h->indices[i] = (int *)calloc(sizeof(int), n_nbrs);
	return (h);
}

double nheap_largest(t_nheap *h, int row)
{
	return (h->distances[row][0]);
}

int nheap_push(t_nheap *h, int row, double val, int i_val)
{
	int i, ic1, ic2, i_swap;
	int size;
	double *dist_arr;
	int *ind_arr;

	size = h->n_nbrs;
	dist_arr = h->distances[row];
	ind_arr = h->indices[row];

	// if distance is already greater than the furthest element, don't push
	if (val > dist_arr[0])
		return (0);

	// insert the values at position 0
	dist_arr[0] = val;
	ind_arr[0] = i_val;

	// descend the heap, swapping values until the max heap criterion is met
	i = 0;
	while (TRUE)
	{
		ic1 = 2 * i + 1;
		ic2 = ic1 + 1;

		if (ic1 >= size)
			break;
		else if (ic2 >= size)
		{
			if (dist_arr[ic1] > val)
				i_swap = ic1;
			else
				break;
		}
		else if (dist_arr[ic1] >= dist_arr[ic2])
		{
			if (val < dist_arr[ic1])
				i_swap = ic1;
			else
				break;
		}
		else
		{
			if (val < dist_arr[ic2])
				i_swap = ic2;
			else
				break;
		}
		dist_arr[i] = dist_arr[i_swap];
		ind_arr[i] = ind_arr[i_swap];
		i = i_swap;
	}

	dist_arr[i] = val;
	ind_arr[i] = i_val;

	return (0);
}

void simultaneous_sort(double *dist, int *idx, int size)
{
	int pivot_idx, i, store_idx;
	double pivot_val;

	if (size <= 1)
		;
	else if (size == 2)
	{
		if (dist[0] > dist[1])
			dual_swap(dist, idx, 0, 1);
	}
	else if (size == 3)
	{
		if (dist[0] > dist[1])
			dual_swap(dist, idx, 0, 1);
		if (dist[1] > dist[2])
		{
			dual_swap(dist, idx, 1, 2);
			if (dist[0] > dist[1])
				dual_swap(dist, idx, 0, 1);
		}
	}
	else
	{
		pivot_idx = size / 2;
		if (dist[0] > dist[size - 1])
			dual_swap(dist, idx, 0, size - 1);
		if (dist[size - 1] > dist[pivot_idx])
		{
			dual_swap(dist, idx, size - 1, pivot_idx);
			if (dist[0] > dist[size - 1])
				dual_swap(dist, idx, 0, size - 1);
		}
		pivot_val = dist[size - 1];

		store_idx = 0;
		for (i = 0; i < size - 1; i++)
		{
			if (dist[i] < pivot_val)
			{
				dual_swap(dist, idx, i, store_idx);
				store_idx++;
			}
		}
		dual_swap(dist, idx, store_idx, size - 1);
		pivot_idx = store_idx;
		if (pivot_idx > 1)
			simultaneous_sort(dist, idx, pivot_idx);
		if (pivot_idx * 2 < size)
			simultaneous_sort(dist + pivot_idx + 1, idx + pivot_idx + 1, size - pivot_idx - 1);
	}
}

void nheap_sort(t_nheap *h)
{
	int row;

	for (row = 0; row < h->n_pts; row++)
		simultaneous_sort(h->distances[row], h->indices[row], h->n_nbrs);
}

t_knn nheap_get_arrays(t_nheap *h)
{
	t_knn output;

	nheap_sort(h);
	output.distances = h->distances;
	output.indices = h->indices;
	output.n_samples = h->n_pts;
	output.n_neighbors = h->n_nbrs;
	free(h);
	return (output);
}

double **copy_double_arr(double **arr, int row, int col)
{
	double **copy;
	int i, j;

	copy = (double **)malloc(sizeof(double *) * row);
	for (i = 0; i < row; i++)
	{
		copy[i] = (double *)malloc(sizeof(double) * col);
		for (j = 0; j < col; j++)
			copy[i][j] = arr[i][j];
	}
	return (copy);
}

void swap_2(int *arr, int i1, int i2)
{
	int tmp;

	tmp = arr[i1];
	arr[i1] = arr[i2];
	arr[i2] = tmp;
}

void btree_zero(t_btree *b)
{
	b->data = NULL;
	b->idx_array = NULL;
	b->node_data = NULL;
	b->node_bounds = NULL;

	b->leaf_size = 40;
	b->n_levels = 0;
	b->n_nodes = 0;
}

int init_node(t_btree *b, int i_node, int idx_start, int idx_end)
{
	int n_features;
	int n_points;
	int i, j;
	double radius;
	double *centroid;

	n_features = b->n_features;
	n_points = idx_end - idx_start;
	centroid = b->node_bounds[0][i_node];

	for (j = 0; j < n_features; j++)
		centroid[j] = 0.0;

	for (i = idx_start; i < idx_end; i++)
		for (j = 0; j < n_features; j++)
			centroid[j] += b->data[b->idx_array[i]][j];

	for (j = 0; j < n_features; j++)
		centroid[j] /= n_points;

	radius = 0.0;
	for (i = idx_start; i < idx_end; i++)
		radius = fmax(radius, manhattan_dist(centroid, b->data[b->idx_array[i]], n_features));

	b->node_data[i_node].radius = radius;
	b->node_data[i_node].idx_start = idx_start;
	b->node_data[i_node].idx_end = idx_end;
	return (0);
}

int find_node_split_dim(double **data, int *node_indices, int n_features, int n_points)
{
	double min_val, max_val, val, spread, max_spread;
	int i, j, j_max;

	j_max = 0;
	max_spread = 0;
	for (j = 0; j < n_features; j++)
	{
		max_val = data[node_indices[0]][j];
		min_val = max_val;
		for (i = 1; i < n_points; i++)
		{
			val = data[node_indices[i]][j];
			max_val = fmax(max_val, val);
			min_val = fmin(min_val, val);
		}
		spread = max_val - min_val;
		if (spread > max_spread)
		{
			max_spread = spread;
			j_max = j;
		}
	}
	return (j_max);
}

int partition_node_indices(double **data, int *node_indices, int split_dim, int split_index,
						   int n_features, int n_points)
{
	(void)n_features;
	int left, right, midindex, i;
	double d1, d2;

	left = 0;
	right = n_points - 1;

	while (TRUE)
	{
		midindex = left;
		for (i = left; i < right; i++)
		{
			d1 = data[node_indices[i]][split_dim];
			d2 = data[node_indices[right]][split_dim];
			if (d1 < d2)
			{
				swap_2(node_indices, i, midindex);
				midindex += 1;
			}
		}
		swap_2(node_indices, midindex, right);
		if (midindex == split_index)
			break;
		else if (midindex < split_index)
			left = midindex + 1;
		else
			right = midindex - 1;
	}
	return (0);
}

void recursive_build(t_btree *b, int i_node, int idx_start, int idx_end)
{
	int imax;
	int n_features;
	int n_points;
	int n_mid;

	n_features = b->n_features;
	n_points = idx_end - idx_start;
	n_mid = n_points / 2;

	// initialize the node data
	init_node(b, i_node, idx_start, idx_end);

	if (2 * i_node + 1 >= b->n_nodes)
	{
		b->node_data[i_node].is_leaf = TRUE;
		if (idx_end - idx_start > 2 * b->leaf_size)
			printf("Memory layout is flawed: not enough nodes allocated");
	}
	else if (idx_end - idx_start < 2)
	{
		printf("Memory layout is flawed: too many nodes allocated");
		b->node_data[i_node].is_leaf = TRUE;
	}
	else
	{
		b->node_data[i_node].is_leaf = FALSE;
		imax = find_node_split_dim(b->data, b->idx_array, n_features, n_points);
		partition_node_indices(b->data, b->idx_array, imax, n_mid, n_features, n_points);
		recursive_build(b, 2 * i_node + 1, idx_start, idx_start + n_mid);
		recursive_build(b, 2 * i_node + 2, idx_start + n_mid, idx_end);
	}
}

t_btree *btree_init(double **data, int n_samples, int n_features, int leaf_size)
{
	t_btree *b;
	int i, j;

	b = (t_btree *)malloc(sizeof(t_btree));
	btree_zero(b);

	b->data = copy_double_arr(data, n_samples, n_features);
	b->leaf_size = leaf_size;

	if (leaf_size < 1)
	{
		printf("leaf_size must be greater than or equal to 1\n");
		exit(-1);
	}

	b->n_samples = n_samples;
	b->n_features = n_features;

	b->n_levels = log2(fmax(1, (b->n_samples - 1) / b->leaf_size)) + 1;
	b->n_nodes = pow(2.0, b->n_levels) - 1;

	b->idx_array = (int *)malloc(sizeof(int) * b->n_samples);
	for (i = 0; i < b->n_samples; i++)
		b->idx_array[i] = i;
	b->node_data = (t_nodedata *)calloc(b->n_nodes, sizeof(t_nodedata));
	b->node_bounds = (double ***)malloc(sizeof(double **));
	b->node_bounds[0] = (double **)malloc(sizeof(double *) * b->n_nodes);
	for (i = 0; i < b->n_nodes; i++)
	{
		b->node_bounds[0][i] = (double *)malloc(sizeof(double) * b->n_features);
		for (j = 0; j < b->n_features; j++)
			b->node_bounds[0][i][j] = 0.0;
	}
	recursive_build(b, 0, 0, b->n_samples);
	return (b);
}

int query_depth_first(t_btree *b, int i_node, double *pt, int i_pt, t_nheap *heap, double dist)
{
	t_nodedata node_info = b->node_data[i_node];
	double dist_pt, dist1, dist2;
	int i, i1, i2;

	// case 1: query point is outside node radius: trim it from the query
	if (dist > nheap_largest(heap, i_pt))
	{
		;
	}
	// case 2: this is a leaf node. Update set of nearby points
	else if (node_info.is_leaf)
	{
		for (i = node_info.idx_start; i < node_info.idx_end; i++)
		{
			dist_pt = manhattan_dist(pt, b->data[b->idx_array[i]], b->n_features);
			if (dist_pt < nheap_largest(heap, i_pt))
				nheap_push(heap, i_pt, dist_pt, b->idx_array[i]);
		}
	}
	// case 3: Node is not a leaf, Recursively query sub-nodes starting with the closest
	else
	{
		i1 = 2 * i_node + 1;
		i2 = i1 + 1;
		dist1 = min_dist(b, i1, pt); // implement min_rdist
		dist2 = min_dist(b, i2, pt);
		if (dist1 <= dist2)
		{
			query_depth_first(b, i1, pt, i_pt, heap, dist1);
			query_depth_first(b, i2, pt, i_pt, heap, dist2);
		}
		else
		{
			query_depth_first(b, i2, pt, i_pt, heap, dist2);
			query_depth_first(b, i1, pt, i_pt, heap, dist1);
		}
	}
	return (0);
}

t_knn btree_query(t_btree *b, double **x, int n_samples, int n_features, int k)
{
	t_nheap *heap;
	double dist;
	int i;
	t_knn output;

	if (n_features != b->n_features)
	{
		printf("query data dimension must match training data dimension.\n");
		exit(-1);
	}
	if (b->n_samples < k)
	{
		printf("k must be less than or equal to the number of training points.\n");
		exit(-1);
	}
	heap = nheap_init(n_samples, k);
	for (i = 0; i < n_samples; i++)
	{
		dist = min_dist(b, 0, x[i]);
		query_depth_first(b, 0, x[i], i, heap, dist);
	}
	output = nheap_get_arrays(heap);
	return (output);
}

void free_2d_double(double **arr, int row)
{
	int i;

	for (i = 0; i < row; i++)
		free(arr[i]);
	free(arr);
}

void free_2d_int(int **arr, int row)
{
	int i;

	for (i = 0; i < row; i++)
		free(arr[i]);
	free(arr);
}

void free_tree(t_btree *tree)
{
	free_2d_double(tree->data, tree->n_samples);
	free(tree->idx_array);
	free(tree->node_data);
	free_2d_double(tree->node_bounds[0], tree->n_nodes);
	free(tree->node_bounds);
	free(tree);
}

void free_knn(t_knn knn, int row)
{
	free_2d_double(knn.distances, row);
	free_2d_int(knn.indices, row);
}
