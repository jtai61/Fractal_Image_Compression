#ifndef BASE_H
#define BASE_H

/* standard library */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <ctype.h>
#include <assert.h>

/* constant define */

#define Baseline		0
#define MassCenter 		1
#define SaupeFisher 	2
#define Saupe 			3
#define Fisher 			4
#define Hurtgen 		5
#define McSaupe 		6
#define Tai				7

#define IDENTITY      	0
#define L_ROTATE90    	1
#define R_ROTATE90    	2
#define ROTATE180     	3
#define R_VERTICAL    	4
#define R_HORIZONTAL  	5
#define F_DIAGONAL    	6
#define S_DIAGONAL    	7

#define MAX_NEIGHBOURS 	1000
#define BUCKETSIZE 		10
#define HEAPMAX 		10000
#define TWOPI 			6.2831853

/* function define */

#define matrix_allocate(matrix, hsize, vsize, TYPE)                           \
	{                                                                         \
		TYPE *imptr;                                                          \
		int _i;                                                               \
		matrix = (TYPE **)malloc((vsize) * sizeof(TYPE *));                   \
		imptr = (TYPE *)malloc((long)(hsize) * (long)(vsize) * sizeof(TYPE)); \
		if (imptr == NULL)                                                    \
			fatal("\nNo memory in matrix allocate.");                         \
		for (_i = 0; _i < vsize; ++_i, imptr += hsize)                        \
			matrix[_i] = imptr;                                               \
	}

#define swap(a, b, TYPE) \
	{                    \
		TYPE _temp;      \
		_temp = b;       \
		b = a;           \
		a = _temp;       \
	}

#define bound(a) ((a) < 0.0 ? 0 : ((a) > 255.0 ? 255 : a))

/* rename data type */

typedef unsigned char PIXEL;

/* structure declaration */

typedef struct kdtree
{
	int cutdim;
	float cutval;
	float *max, *min;
	struct kdtree *left, *right;
	int *idx, num;
} kdtree;

typedef struct heapT
{
	float dist;
	kdtree *t;
} heapT;

struct code_book
{
	short ptr_x;
	short ptr_y;
	short isom;
	double sum;
	double sum2;
};

struct c
{
	short ptr_x;
	short ptr_y;
	double sum;
	double sum2;
	short iso;
	struct c *next;
};

struct t_node
{
	double rx;
	double ry;
	double rrx;
	double rry;
	short size;
	double dx;
	double dy;
	short sym_op;
	double alfa, beta;
	struct t_node *next;
};

/* global variable */

#ifndef EXTERN
#define EXTERN extern
#define INIT(x)
extern int ordering[][6];
extern int mapping[][8];
#else
#define INIT(x) x

int ordering[][6] = {{3, 2, 1, 0, 3, 0},
					 {2, 3, 1, 0, 5, 1},
					 {2, 1, 3, 0, 5, 2},
					 {3, 1, 2, 0, 7, 0},
					 {1, 3, 2, 0, 1, 1},
					 {1, 2, 3, 0, 1, 2},
					 {3, 2, 0, 1, 3, 1},
					 {2, 3, 0, 1, 5, 0},
					 {2, 1, 0, 3, 2, 2},
					 {3, 1, 0, 2, 7, 1},
					 {1, 3, 0, 2, 1, 0},
					 {1, 2, 0, 3, 4, 2},
					 {3, 0, 2, 1, 3, 2},
					 {2, 0, 3, 1, 2, 0},
					 {2, 0, 1, 3, 2, 1},
					 {3, 0, 1, 2, 7, 2},
					 {1, 0, 2, 3, 4, 1},
					 {1, 0, 3, 2, 4, 0},
					 {0, 3, 2, 1, 6, 2},
					 {0, 2, 3, 1, 6, 1},
					 {0, 2, 1, 3, 6, 0},
					 {0, 3, 1, 2, 0, 2},
					 {0, 1, 3, 2, 0, 1},
					 {0, 1, 2, 3, 0, 0}};

int mapping[][8] = {{0, 1, 2, 3, 4, 5, 6, 7},
					{2, 0, 3, 1, 7, 6, 4, 5},
					{1, 3, 0, 2, 6, 7, 5, 4},
					{3, 2, 1, 0, 5, 4, 7, 6},
					{4, 7, 6, 5, 0, 3, 2, 1},
					{5, 6, 7, 4, 3, 0, 1, 2},
					{6, 4, 5, 7, 1, 2, 0, 3},
					{7, 5, 4, 6, 2, 1, 3, 0}};
#endif

EXTERN struct t_node *trans, fractal_code;

EXTERN int min;
EXTERN int max;
EXTERN int lev;

EXTERN FILE *input;
EXTERN FILE *output;
EXTERN int iterations INIT(= 8);
EXTERN int postproc INIT(= 0);
EXTERN int quality INIT(= 0);
EXTERN int piramidal INIT(= 1);
EXTERN int raw_format INIT(= 0);
EXTERN int display INIT(= 0);
EXTERN double zoom INIT(= 1.0);

EXTERN PIXEL **image1;

EXTERN struct c ***class_polar[8];
EXTERN struct c *class_fisher[8][3][24];
EXTERN struct c *class_hurtgen[8][16][24];

EXTERN kdtree ***class_polar_saupe[8];
EXTERN int **item_per_class[8];
EXTERN struct code_book ***c_book[8];
EXTERN float ****f_vect[8];

EXTERN struct code_book *codebook[8];
EXTERN struct code_book *codebook_v2[8][3];
EXTERN float **f_vectors[8];
EXTERN float **f_vectors_v2[8][3];
EXTERN kdtree *kd_tree[8];
EXTERN kdtree *kd_tree_v2[8][3];
EXTERN int feat_vect_dim[8];
EXTERN int average_factor[8];
EXTERN int clas_count[8][3];

EXTERN clock_t start_clock;
EXTERN clock_t end_clock;

EXTERN int n_features INIT(= 16);
EXTERN int shrunk_factor_mc INIT(= 1);
EXTERN int shrunk_factor_saupe INIT(= 0);

EXTERN int N_BITALFA INIT(= 4);
EXTERN int N_BITBETA INIT(= 7);
EXTERN double MAX_ALFA INIT(= 1.0);
EXTERN int min_size INIT(= 4);
EXTERN int max_size INIT(= 16);
EXTERN int image_width INIT(= 512);
EXTERN int image_height INIT(= 512);
EXTERN int virtual_size INIT(= 512);
EXTERN int SHIFT INIT(= 4);

EXTERN double T_ENT INIT(= 8.0);
EXTERN double T_VAR INIT(= 1000000.0);
EXTERN double T_RMS INIT(= 8.0);
EXTERN double adapt INIT(= 1.0);

EXTERN double eps INIT(= 2.0);
EXTERN int matches INIT(= 50);
EXTERN int n_p_class INIT(= 50);
EXTERN int livelli INIT(= 3);
EXTERN int zero_threshold INIT(= 0);
EXTERN int best INIT(= 0);
EXTERN int full_first_class INIT(= 0);
EXTERN int full_second_class INIT(= 0);
EXTERN int qtree INIT(= 0);
EXTERN FILE *fp;

EXTERN int bits_per_coordinate_w;
EXTERN int bits_per_coordinate_h;

EXTERN int zero_alfa_transform INIT(= 0);
EXTERN long comparisons INIT(= 0);
EXTERN int transforms INIT(= 0);
EXTERN int zeroalfa;

EXTERN char filein[100];
EXTERN char fileout[100];

EXTERN PIXEL **image;
EXTERN PIXEL **qtt;
EXTERN double **contract;
EXTERN double **range;
EXTERN double **range_tmp;
EXTERN double **flip_range;
EXTERN double (*Coding)(int, int, int, int *, int *, int *, int *, int *);
EXTERN void (*Indexing)(int, int);
EXTERN int method INIT(= MassCenter);

/* function prototype */

double BaselineCoding(int, int, int, int *, int *, int *, int *, int *);
double HurtgenCoding(int, int, int, int *, int *, int *, int *, int *);
double SaupeCoding(int, int, int, int *, int *, int *, int *, int *);
double FisherCoding(int, int, int, int *, int *, int *, int *, int *);
double Mc_SaupeCoding(int, int, int, int *, int *, int *, int *, int *);
double MassCenterCoding(int, int, int, int *, int *, int *, int *, int *);
double Saupe_FisherCoding(int, int, int, int *, int *, int *, int *, int *);
double TaiCoding(int, int, int, int *, int *, int *, int *, int *);
double entropy(int, int, int, int);
double variance(int, int, int, int);
double variance_2(int, double **, int, int);
int hurtgen_class(int, double **);
int pack(int, long, FILE *);
int variance_class(int, double **);
void ComputeAverageFactorMc();
void ComputeFeatVectDimSaupe();
void ComputeMcVectors(double **, double **, int, int, double *);
void ComputeSaupeVectors(double **, int, int, float *);
void BaselineIndexing(int, int);
void FisherIndexing(int, int);
void HurtgenIndexing(int, int);
void MassCenterIndexing(int, int);
void SaupeIndexing(int, int);
void Mc_SaupeIndexing(int, int);
void Saupe_FisherIndexing(int, int);
void TaiIndexing(int, int);
void contraction(double **, PIXEL **, int, int);
void fatal(const char *);
void flips(int, double **, double **, int);
void ComputeMc(double **, int, double *, double *, int);
void newclass(int, double **, int *, int *);
void newclass_v2(int, double **, int *, int *);
void getopt_enc(int, char **);
void getopt_dec(int, char **);
void quadtree(int, int, int, double, double, double);
void quadtree_v2(int, int, int, double, double, double);
void readimage_raw(char *);
void readimage_pgm(char *, int *, int *);
void help_enc();
kdtree *kdtree_build(float **, int, int);
kdtree *kdtree_build_v2(float **, int, int);
void kdtree_free(kdtree *);
int kdtree_search(float *, float **, int, kdtree *, float, int, int *);
int kdtree_search_v2(float *, float **, int, kdtree *, float, int, int *);
long unpack(int, FILE *);
void read_transformations(int, int, int);
void writeimage_pgm(const char *, PIXEL **, int, int);
void writeimage_raw(const char *, PIXEL **, int, int);
void writeimage_pipe(FILE *, PIXEL **, int, int);
void smooth_image();
void zooming(double);
void help_dec();
void iterative_decoding(int, int, double);
void piramidal_decoding(int);
double calc_PSNR(char *, char *);
double calc_SSIM(char *, char *);

#endif