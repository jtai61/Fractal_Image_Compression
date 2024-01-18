#include <stdio.h>
#include <math.h>
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
