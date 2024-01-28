#include <stdio.h>
#include "nns.h"
#include "def.h"

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

typedef unsigned char PIXEL;
typedef unsigned int LBP;

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

EXTERN struct t_node *trans, fractal_code;

EXTERN int min;
EXTERN int max;
EXTERN int lev;

EXTERN FILE *input;
EXTERN FILE *output;
EXTERN int iterations INIT(= 10);
EXTERN int postproc INIT(= 0);
EXTERN int quality INIT(= 0);
EXTERN int piramidal INIT(= 1);
EXTERN int raw_format INIT(= 0);
EXTERN int display INIT(= 0);
EXTERN double zoom INIT(= 1.0);

EXTERN PIXEL **image1;

EXTERN struct c ***class_polar[8];
EXTERN struct c *class_fisher[8][3][24];
EXTERN struct c *class_nandi[8][3][24];
EXTERN struct c *class_tai[8][3][24];
EXTERN struct c *class_hurtgen[8][16][24];

EXTERN kdtree ***class_polar_saupe[8];
EXTERN int **item_per_class[8];
EXTERN struct code_book ***c_book[8];
EXTERN float ****f_vect[8];

EXTERN struct code_book *codebook[8];
EXTERN float **f_vectors[8];
EXTERN kdtree *kd_tree[8];
EXTERN int feat_vect_dim[8];
EXTERN int average_factor[8];
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

EXTERN char filein[50];
EXTERN char fileout[50];

EXTERN PIXEL **image;
EXTERN PIXEL **qtt;
EXTERN double **contract;
EXTERN double **range;
EXTERN double **range_tmp;
EXTERN double **flip_range;
EXTERN double (*Coding)();
EXTERN void (*Indexing)();
EXTERN int method INIT(= MassCenter);
