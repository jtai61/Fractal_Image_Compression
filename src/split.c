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

    left_up_pixel = block[xFloor][yFloor];
    right_up_pixel = block[xFloor][yFloor + 1];
    left_down_pixel = block[xFloor + 1][yFloor];
    right_down_pixel = block[xFloor + 1][yFloor + 1];

    // interpolation weights
    float alpha, beta;

    alpha = target_y - yFloor;
    beta = target_x - xFloor;

    return (1 - alpha) * (1 - beta) * left_up_pixel + alpha * (1 - beta) * right_up_pixel + (1 - alpha) * beta * left_down_pixel + alpha * beta * right_down_pixel;
}

int HammingDistance(LBP num1, LBP num2)
{
    LBP num_xor = num1 ^ num2;
    int distance = 0;

    while (num_xor)
    {
        ++distance;
        num_xor &= (num_xor - 1);
    }

    return distance;
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
