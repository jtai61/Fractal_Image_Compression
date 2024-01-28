#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "def.h"
#include "globals.h"

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

int HammingDistance(LBP lbp1, LBP lbp2)
{
    LBP lbp1_lbp2_xor = lbp1 ^ lbp2;
    int distance = 0;

    /* Brian Kernighan's algorithm */
    while (lbp1_lbp2_xor)
    {
        ++distance;
        lbp1_lbp2_xor &= (lbp1_lbp2_xor - 1);
    }

    return distance;
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