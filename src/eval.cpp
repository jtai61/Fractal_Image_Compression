#define EXTERN
#include "base.hpp"

int main(int argc, char **argv)
{
    printf("\nPSNR : %.2f dB\n\n", calc_PSNR(argv[1], argv[2]));
    printf("\nSSIM : %.4f\n", calc_SSIM(argv[1], argv[2]));


    return 0;
}

double calc_PSNR(char *img1_path, char *img2_path)
{
    FILE *in1, *in2;

    if ((in1 = fopen(img1_path, "rb")) == NULL || (in2 = fopen(img2_path, "rb")) == NULL)
        fatal("Can't open input file !\n");

    PIXEL **img1, **img2;

    matrix_allocate(img1, image_width, image_height, PIXEL);
    matrix_allocate(img2, image_width, image_height, PIXEL);

    printf("Reading %s (%dx%d) ...\n", img1_path, image_width, image_height);
    fflush(stdout);
    printf("Reading %s (%dx%d) ...\n", img2_path, image_width, image_height);
    fflush(stdout);

    int i, j;

    i = fread(img1[0], sizeof(PIXEL), image_width * image_height, in1);
    j = fread(img2[0], sizeof(PIXEL), image_width * image_height, in2);

    if (i < image_width * image_height || j < image_width * image_height)
        fatal("Not enough input data in the input file !\n");

    fclose(in1);
    fclose(in2);

    printf("\nStart calculate PSNR ...\n");

    int max_pixel_value = 255;
    double sum = 0.0, mse, psnr;

    for (i = 0; i < image_height; i++)
        for (j = 0; j < image_width; j++)
            sum += pow(img1[i][j] - img2[i][j], 2);

    mse = sum / (double)(image_width * image_height);

    if (mse == 0.0)
        psnr = INFINITY;
    else
        psnr = 10.0 * log10(pow(max_pixel_value, 2) / mse);

    free(img1[0]);
    free(img2[0]);

    return psnr;
}

double calc_SSIM(char *img1_path, char *img2_path)
{
    FILE *in1, *in2;

    if ((in1 = fopen(img1_path, "rb")) == NULL || (in2 = fopen(img2_path, "rb")) == NULL)
        fatal("Can't open input file !\n");

    PIXEL **img1, **img2;

    matrix_allocate(img1, image_width, image_height, PIXEL);
    matrix_allocate(img2, image_width, image_height, PIXEL);

    printf("Reading %s (%dx%d) ...\n", img1_path, image_width, image_height);
    fflush(stdout);
    printf("Reading %s (%dx%d) ...\n", img2_path, image_width, image_height);
    fflush(stdout);

    int i, j;

    i = fread(img1[0], sizeof(PIXEL), image_width * image_height, in1);
    j = fread(img2[0], sizeof(PIXEL), image_width * image_height, in2);

    if (i < image_width * image_height || j < image_width * image_height)
        fatal("Not enough input data in the input file !\n");

    fclose(in1);
    fclose(in2);

    printf("\nStart calculate SSIM ...\n");

    /* calculate mean */
    double img1_mean = 0.0, img2_mean = 0.0;

    for (i = 0; i < image_height; i++)
        for (j = 0; j < image_width; j++)
            img1_mean += img1[i][j];

    img1_mean /= (double)(image_height * image_width);

    for (i = 0; i < image_height; i++)
        for (j = 0; j < image_width; j++)
            img2_mean += img2[i][j];

    img2_mean /= (double)(image_height * image_width);

    /* calculate standard deviation */
    double img1_std = 0.0, img2_std = 0.0;

    for (i = 0; i < image_height; i++)
        for (j = 0; j < image_width; j++)
            img1_std += pow(img1[i][j] - img1_mean, 2);

    img1_std = sqrt(img1_std / (double)(image_height * image_width));

    for (i = 0; i < image_height; i++)
        for (j = 0; j < image_width; j++)
            img2_std += pow(img2[i][j] - img2_mean, 2);

    img2_std = sqrt(img2_std / (double)(image_height * image_width));

    /* calculate covariance */
    double img1_img2_cov = 0.0;

    for (i = 0; i < image_height; i++)
        for (j = 0; j < image_width; j++)
            img1_img2_cov += (img1[i][j] - img1_mean) * (img2[i][j] - img2_mean);

    img1_img2_cov /= (double)(image_height * image_width);

    /* calculate SSIM */
    const double K1 = 0.01, K2 = 0.03, L = 255.0;
    const double C1 = pow(K1 * L, 2), C2 = pow(K2 * L, 2), C3 = C2 / 2.0;

    double l, c, s, ssim;

    l = (2.0 * img1_mean * img2_mean + C1) / (pow(img1_mean, 2) + pow(img2_mean, 2) + C1);
    c = (2.0 * img1_std * img2_std + C2) / (pow(img1_std, 2) + pow(img2_std, 2) + C2);
    s = (img1_img2_cov + C3) / (img1_std * img2_std + C3);

    ssim = l * c * s;

    free(img1[0]);
    free(img2[0]);

    return ssim;
}