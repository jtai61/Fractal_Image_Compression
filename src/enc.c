#define EXTERN
#include "base.h"

int main(int argc, char **argv)
{
    start_clock = clock();

    int i, h, k, max;
    char *p, *pp;
    int int_max_alfa;

    getopt_enc(argc, argv);

    pp = (char *)"";
    p = filein;
    while (*p)
    { /* Find the last dot in the file name */
        if (*p == '.')
            pp = p + 1;
        p++;
    }

    if (!strcmp(pp, "pgm"))
        readimage_pgm(filein, &image_width, &image_height);
    else
        readimage_raw(filein);

    max = image_height;
    if (image_width > image_height)
        max = image_width;

    virtual_size = 1 << (int)ceil(log((double)max) / log(2.0));

    matrix_allocate(contract, 1 + image_width / 2, 1 + image_height / 2, double);
    matrix_allocate(qtt, virtual_size, virtual_size, PIXEL);
    matrix_allocate(range_tmp, 64, 64, double);
    matrix_allocate(flip_range, 64, 64, double);
    matrix_allocate(range, 64, 64, double);

    switch (method)
    {
    case Baseline:
        Indexing = BaselineIndexing;
        Coding = BaselineCoding;
        printf("\n Speed-up method: Baseline\n\n");
        break;

    case MassCenter:
        ComputeAverageFactorMc();
        for (k = 0; k < 8; k++)
        {
            matrix_allocate(class_polar[k], n_p_class, n_p_class, struct c *);
        }

        for (k = 0; k < 8; k++)
            for (h = 0; h < n_p_class; h++)
                for (i = 0; i < n_p_class; i++)
                    class_polar[k][h][i] = NULL;

        Indexing = MassCenterIndexing;
        Coding = MassCenterCoding;
        printf("\n Speed-up method: MassCenter\n\n");
        break;

    case SaupeFisher:
        ComputeFeatVectDimSaupe();
        Indexing = Saupe_FisherIndexing;
        Coding = Saupe_FisherCoding;
        printf("\n Speed-up method: Saupe-Fisher\n\n");
        break;

    case Saupe:
        ComputeFeatVectDimSaupe();
        Indexing = SaupeIndexing;
        Coding = SaupeCoding;
        printf("\n Speed-up method: Saupe\n\n");
        break;

    case Fisher:
        for (k = 0; k < 8; k++)
            for (h = 0; h < 3; h++)
                for (i = 0; i < 24; i++)
                    class_fisher[k][h][i] = NULL;

        Indexing = FisherIndexing;
        Coding = FisherCoding;
        printf("\n Speed-up method: Fisher\n\n");
        break;

    case Hurtgen:
        for (k = 0; k < 8; k++)
            for (h = 0; h < 16; h++)
                for (i = 0; i < 24; i++)
                    class_hurtgen[k][h][i] = NULL;
        Indexing = HurtgenIndexing;
        Coding = HurtgenCoding;
        printf("\n Speed-up method: Hurtgen\n\n");
        break;

    case McSaupe:
        ComputeAverageFactorMc();
        ComputeFeatVectDimSaupe();

        for (k = 0; k < 8; k++)
        {
            matrix_allocate(class_polar[k], n_p_class, n_p_class, struct c *);
            matrix_allocate(class_polar_saupe[k], n_p_class, n_p_class, kdtree *);
            matrix_allocate(item_per_class[k], n_p_class, n_p_class, int);
            matrix_allocate(c_book[k], n_p_class, n_p_class, struct code_book *);
            matrix_allocate(f_vect[k], n_p_class, n_p_class, float **);
        }

        for (k = 0; k < 8; k++)
            for (h = 0; h < n_p_class; h++)
                for (i = 0; i < n_p_class; i++)
                {
                    class_polar[k][h][i] = NULL;
                    class_polar_saupe[k][h][i] = NULL;
                    item_per_class[k][h][i] = 0;
                }
        Indexing = Mc_SaupeIndexing;
        Coding = Mc_SaupeCoding;
        printf("\n Speed-up method: Mc-Saupe\n\n");
        break;

    case Tai:
        ComputeFeatVectDimSaupe();
        Indexing = TaiIndexing;
        Coding = TaiCoding;
        printf("\n Speed-up method: Tai\n\n");
        break;
    }

    contraction(contract, image, 0, 0);

    for (i = (int)rint(log((double)2) / log(2.0)); i <= (int)rint(log((double)max_size) / log(2.0)); i++)
        Indexing((int)rint(pow(2.0, (double)i)), i);

    bits_per_coordinate_w = ceil(log(image_width / SHIFT) / log(2.0));
    bits_per_coordinate_h = ceil(log(image_height / SHIFT) / log(2.0));

    for (k = 0; k < virtual_size; k++)
        for (h = 0; h < virtual_size; h++)
            qtt[k][h] = 255;

    for (k = 0; k < virtual_size; k++)
    {
        qtt[k][0] = 0;
        qtt[k][virtual_size - 1] = 0;
    }

    for (k = 0; k < virtual_size; k++)
    {
        qtt[0][k] = 0;
        qtt[virtual_size - 1][k] = 0;
    }

    int_max_alfa = 0.5 + (MAX_ALFA) / (8.0) * (1 << 8);
    if (int_max_alfa < 0)
        int_max_alfa = 0;
    if (int_max_alfa >= (1 << 8))
        int_max_alfa = (1 << 8) - 1;

    MAX_ALFA = (double)int_max_alfa / (double)(1 << 8) * (8.0);
    zeroalfa = 0;

    if ((fp = fopen(fileout, "w")) == NULL)
        fatal("\nCan't open output file");

    /* Header of output file */
    pack(4, (long)N_BITALFA, fp);
    pack(4, (long)N_BITBETA, fp);
    pack(7, (long)min_size, fp);
    pack(7, (long)max_size, fp);
    pack(6, (long)SHIFT, fp);
    pack(12, (long)image_width, fp);
    pack(12, (long)image_height, fp);
    pack(8, (long)int_max_alfa, fp);

    printf("\n Image Entropy      : %.2f\n", entropy(image_width, image_height, 0, 0));
    printf(" Image Variance     : %.2f\n", variance(image_width, image_height, 0, 0));
    printf(" Entropy threshold  : %.2f\n", T_ENT);
    printf(" Variance threshold : %.2f\n", T_VAR);
    printf(" Rms threshold      : %.2f\n\n", T_RMS);

    if (method == Tai)
    {
        quadtree_v2(0, 0, virtual_size, T_ENT, T_RMS, T_VAR);
    }
    else
    {
        quadtree(0, 0, virtual_size, T_ENT, T_RMS, T_VAR);
    }

    pack(-1, (long)0, fp);
    i = pack(-2, (long)0, fp);
    fclose(fp);

    printf("\n\n Zero_alfa_transformations   : %d\n", zero_alfa_transform);
    printf(" Number of transformations   : %d\n", transforms);
    printf(" Number of comparisons       : %ld\n", comparisons);
    printf(" Comparisons/Transformations : %.2f\n", (double)comparisons / transforms);
    printf(" %d bytes written in %s\n", i, fileout);

    if (qtree)
        writeimage_pgm("quadtree.pgm", qtt, image_width, image_height);

    free(image[0]);
    free(contract[0]);
    free(flip_range[0]);
    free(range_tmp[0]);
    free(range[0]);

    end_clock = clock();

    printf("\n encode time : %.4f sec\n\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);


    return (0);
}
