#define IDENTITY      0      /*                          */
#define L_ROTATE90    1      /*                          */
#define R_ROTATE90    2      /*                          */
#define ROTATE180     3      /*       Isometries         */
#define R_VERTICAL    4      /*                          */
#define R_HORIZONTAL  5      /*                          */
#define F_DIAGONAL    6      /*                          */
#define S_DIAGONAL    7      /*                          */

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

#define MAX_NEIGHBOURS 1000

#define MassCenter     0
#define SaupeFisher    1
#define Saupe          2
#define Fisher         3
#define Hurtgen        4
#define McSaupe        5 

#define TWOPI 6.2831853 

