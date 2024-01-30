#include "base.h"

double entropy(int width, int height, int atx, int aty)
{
    int i, j;
    double frequency[256];
    double element[256];
    double entrop = 0.0;
    double number;

    number = (width * height);

    for (i = 0; i < 256; i++)
        element[i] = 0.0;

    for (i = atx; i < atx + height; i++)
        for (j = aty; j < aty + width; j++)
        {
            element[image[i][j]]++;
        }

    for (i = 0; i < 256; i++)
        frequency[i] = element[i] / number;

    for (i = 0; i < 256; i++)
        if (frequency[i] > 0)
            entrop += frequency[i] * log((double)(1.0 / frequency[i])) / log(2.0);

    return entrop;
}

double variance(int width, int height, int atx, int aty)
{
    int i, j;
    double sum2 = 0.0;
    double sum = 0.0;
    double mean2;
    double mean;

    for (i = atx; i < atx + height; i++)
        for (j = aty; j < aty + width; j++)
        {
            sum2 += image[i][j] * image[i][j];
            sum += image[i][j];
        }

    mean2 = sum2 / (width * height);
    mean = sum / (width * height);

    return mean2 - mean * mean;
}

double variance_2(int size, double **block, int atx, int aty)
{
    int i, j;
    double sum2 = 0.0;
    double sum = 0.0;
    double mean2;
    double mean;

    for (i = atx; i < atx + size; i++)
        for (j = aty; j < aty + size; j++)
        {
            sum2 += block[i][j] * block[i][j];
            sum += block[i][j];
        }

    mean2 = sum2 / (size * size);
    mean = sum / (size * size);

    return mean2 - mean * mean;
}

double mean(int size, double **block, int atx, int aty)
{
    int i, j;
    double sum = 0.0;

    for (i = atx; i < atx + size; i++)
        for (j = aty; j < aty + size; j++)
            sum += block[i][j];

    return sum / (double)(size * size);
}

double BilinearInterpolation(double **block, float target_x, float target_y)
{
    int xFloor = (int)target_x, yFloor = (int)target_y;

    // four nearest pixels
    double left_up_pixel, right_up_pixel, left_down_pixel, right_down_pixel;

    left_up_pixel = block[yFloor][xFloor];
    right_up_pixel = block[yFloor][xFloor + 1];
    left_down_pixel = block[yFloor + 1][xFloor];
    right_down_pixel = block[yFloor + 1][xFloor + 1];

    // interpolation weights
    float alpha, beta;

    alpha = target_x - xFloor;
    beta = target_y - yFloor;

    return (1 - alpha) * (1 - beta) * left_up_pixel + alpha * (1 - beta) * right_up_pixel + (1 - alpha) * beta * left_down_pixel + alpha * beta * right_down_pixel;
}

int ELBP_HammingDistance(LBP ci_1, LBP ci_2, LBP ni_1, LBP ni_2, LBP rd_1, LBP rd_2)
{
    LBP ci_xor = ci_1 ^ ci_2, ni_xor = ni_1 ^ ni_2, rd_xor = rd_1 ^ rd_2;
    int ci_hd = 0, ni_hd = 0, rd_hd = 0;

    /* calculate two ELBP_CI hamming distance */
    while (ci_xor)
    {
        ++ci_hd;
        ci_xor &= (ci_xor - 1);
    }

    /* calculate two ELBP_NI hamming distance */
    while (ni_xor)
    {
        ++ni_hd;
        ni_xor &= (ni_xor - 1);
    }

    /* calculate two ELBP_RD hamming distance */
    while (rd_xor)
    {
        ++rd_hd;
        rd_xor &= (rd_xor - 1);
    }

    return ci_hd + ni_hd + rd_hd;
}

void countingSort(int *data, int size)
{
    if (size <= 1) /* sorted */
        return;
    else
    {
        int i;

        /* find max */
        int max = data[0];

        for (i = 1; i < size; i++)
        {
            if (data[i] > max)
                max = data[i];
        }

        /* create counting array */
        int count[max + 1];

        for (i = 0; i <= max; i++)
            count[i] = 0;

        /* count the number of each element */
        for (i = 0; i < size; i++)
            count[data[i]]++;

        /* reconstruct sorted array */
        int index = 0;

        for (i = 0; i <= max; i++)
        {
            while (count[i] > 0)
            {
                data[index++] = i;
                count[i]--;
            }
        }
    }
}
