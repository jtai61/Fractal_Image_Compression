#define EXTERN
#include "base.h"

int main(int argc, char **argv)
{
	start_clock = clock();
	
	int int_max_alfa, step;
	int pipe_disp[2], pid;
	FILE *output;

	getopt_dec(argc, argv);

	if ((input = fopen(filein, "r")) == NULL)
		fatal("\n Can't open input file");

	unpack(-2, input); /* Initialize unpack */

	N_BITALFA = (int)unpack(4, input);
	N_BITBETA = (int)unpack(4, input);
	min_size = (int)unpack(7, input);
	max_size = (int)unpack(7, input);
	SHIFT = (int)unpack(6, input);
	image_width = (int)unpack(12, input);
	image_height = (int)unpack(12, input);
	int_max_alfa = (int)unpack(8, input);

	bits_per_coordinate_w = ceil(log(image_width / SHIFT) / log(2.0));
	bits_per_coordinate_h = ceil(log(image_height / SHIFT) / log(2.0));

	zeroalfa = 0;
	MAX_ALFA = (double)int_max_alfa / (double)(1 << 8) * (8.0);

	max = image_height;
	min = image_width;
	if (image_width > image_height)
	{
		min = image_height;
		max = image_width;
	}

	virtual_size = 1 << (int)ceil(log((double)max) / log(2.0));

	trans = &fractal_code;
	printf("\n Reading %s ... ", filein);
	fflush(stdout);
	read_transformations(0, 0, virtual_size);
	printf("done\n");
	fflush(stdout);

	printf(" Original image size: %dx%d\n", image_width, image_height);

	image_width = (int)rint((zoom * image_width));
	image_height = (int)rint((zoom * image_height));

	if (zoom != 1.0)
	{
		printf(" Zooming image to   : %dx%d\n", image_width, image_height);
		fflush(stdout);
	}

	matrix_allocate(image, 2 + image_width, 2 + image_height, PIXEL);
	matrix_allocate(image1, 2 + image_width, 2 + image_height, PIXEL);

	if (piramidal)
	{
		min *= zoom;
		step = SHIFT * floor(zoom);
		if (step == 0)
			step = 1;
		lev = 0;
		while (1)
		{
			if (min < 200 || (step & 1))
				break;
			min >>= 1;
			step >>= 1;
			lev++;
		}
		printf("\n %d level piramid\n", lev);

		iterative_decoding(lev, iterations, zoom); /* Decode at low resolution */
		piramidal_decoding(lev);				   /* Increase resolution      */
		if (quality)
			iterative_decoding(0, 2, 1.0);
	}
	else
		iterative_decoding(0, iterations, zoom);

	if (postproc)
		smooth_image();

	if (display)
	{
		if (pipe(pipe_disp))
			fatal("\n Pipe creation error");

		if ((pid = fork()) < 0)
			fatal("\n Fork error");

		if (pid == 0)
		{ /* Child */
			dup2(pipe_disp[0], fileno(stdin));
			close(pipe_disp[0]);
			close(pipe_disp[1]);
			execlp("xv", "xv", "-", (char *)0);
			fatal("\n Exec error xv not started");
		}

		printf("\n");
		close(pipe_disp[0]);
		output = fdopen(pipe_disp[1], "w");
		writeimage_pipe(output, image, image_width, image_height);
	}
	else if (raw_format)
		writeimage_raw(fileout, image, image_width, image_height);
	else
		writeimage_pgm(fileout, image, image_width, image_height);

	free(image[0]);
	free(image1[0]);

	end_clock = clock();

	printf("\n PSNR : %.2f", calc_PSNR("img/uncompressed/Lena_512.raw", "img/decompressed/Lena_512_dec.raw"));
	printf("\n SSIM : %.4f\n", calc_SSIM("img/uncompressed/Lena_512.raw", "img/decompressed/Lena_512_dec.raw"));

	printf("\n decode time : %.4f sec\n\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);

	return 0;
}

void zooming(double scalefactor)
{
	trans = &fractal_code;
	while (trans->next != NULL)
	{
		trans = trans->next;

		trans->rrx *= scalefactor;
		trans->rry *= scalefactor;
		trans->rx *= scalefactor;
		trans->ry *= scalefactor;
		trans->dx *= scalefactor;
		trans->dy *= scalefactor;
	}
}

void read_transformations(int atx, int aty, int size)
{
	int qalfa, qbeta;
	double alfa, beta;

	if (atx >= image_height || aty >= image_width)
		return;

	if (size > max_size || atx + size > image_height || aty + size > image_width)
	{
		read_transformations(atx, aty, size / 2);
		read_transformations(atx + size / 2, aty, size / 2);
		read_transformations(atx, aty + size / 2, size / 2);
		read_transformations(atx + size / 2, aty + size / 2, size / 2);
		return;
	}

	if (size > min_size && unpack(1, input))
	{
		/* A 1 means we subdivided.. so quadtree */
		read_transformations(atx, aty, size / 2);
		read_transformations(atx + size / 2, aty, size / 2);
		read_transformations(atx, aty + size / 2, size / 2);
		read_transformations(atx + size / 2, aty + size / 2, size / 2);
	}
	else
	{
		/* Read the trasformation */
		trans->next = (struct t_node *)malloc(sizeof(struct t_node));
		trans = trans->next;
		trans->next = NULL;
		qalfa = (int)unpack(N_BITALFA, input);
		qbeta = (int)unpack(N_BITBETA, input);

		/* Compute alfa from the quantized value */
		alfa = (double)qalfa / (double)(1 << N_BITALFA) * (MAX_ALFA);

		/* Compute beta from the quantized value */
		beta = (double)qbeta / (double)((1 << N_BITBETA) - 1) * ((1.0 + fabs(alfa)) * 255);
		if (alfa > 0.0)
			beta -= alfa * 255;

		trans->alfa = alfa;
		trans->beta = beta;
		if (qalfa != zeroalfa)
		{
			trans->sym_op = (int)unpack(3, input);
			trans->dx = SHIFT * (int)unpack(bits_per_coordinate_h, input);
			trans->dy = SHIFT * (int)unpack(bits_per_coordinate_w, input);
		}
		else
		{
			trans->sym_op = 0;
			trans->dx = 0;
			trans->dy = 0;
		}
		trans->rx = atx;
		trans->ry = aty;
		trans->size = size;
		trans->rrx = atx + size;
		trans->rry = aty + size;
	}
}

void iterative_decoding(int level, int n_iter, double zoo)
{
	int rx, ry, rrx, rry, dx, dy;
	int i, j, ii, jj, s;
	PIXEL **imag, **imag1;
	double pixel;
	double z_factor;
	int width, height;
	static int first_time = 0;

	printf("\n");
	z_factor = zoo / (double)(1 << level);
	zooming(z_factor);
	width = (int)rint(image_width * z_factor / zoo);
	height = (int)rint(image_height * z_factor / zoo);

	if (first_time++ == 0)
		for (i = 0; i < height; i++)
			for (j = 0; j < width; j++)
				image[i][j] = 0;

	for (s = 0; s < n_iter; s++)
	{
		imag = image;
		imag1 = image1;

		if (level > 0)
			printf(" Decoding at low resolution (%dx%d) %d\r", width, height, s);
		else
			printf(" Iterative decoding... %d\r", s);

		fflush(stdout);
		trans = &fractal_code;
		while (trans->next != NULL)
		{
			trans = trans->next;
			rx = (int)rint(trans->rx);
			ry = (int)rint(trans->ry);
			dx = (int)rint(trans->dx);
			dy = (int)rint(trans->dy);
			rrx = (int)rint(trans->rrx);
			rry = (int)rint(trans->rry);

			switch (trans->sym_op)
			{
			case IDENTITY:
				for (i = rx, ii = dx; i < rrx; i++, ii += 2)
					for (j = ry, jj = dy; j < rry; j++, jj += 2)
					{
						pixel = (double)(imag[ii][jj] + imag[ii + 1][jj] + imag[ii][jj + 1] + imag[ii + 1][jj + 1]) / 4.0;
						imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
					}
				break;
			case R_ROTATE90:
				for (j = rry - 1, ii = dx; j >= ry; j--, ii += 2)
					for (i = rx, jj = dy; i < rrx; i++, jj += 2)
					{
						pixel = (double)(imag[ii][jj] + imag[ii + 1][jj] + imag[ii][jj + 1] + imag[ii + 1][jj + 1]) / 4.0;
						imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
					}
				break;
			case L_ROTATE90:
				for (j = ry, ii = dx; j < rry; j++, ii += 2)
					for (i = rrx - 1, jj = dy; i >= rx; i--, jj += 2)
					{
						pixel = (double)(imag[ii][jj] + imag[ii + 1][jj] + imag[ii][jj + 1] + imag[ii + 1][jj + 1]) / 4.0;
						imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
					}
				break;
			case ROTATE180:
				for (i = rrx - 1, ii = dx; i >= rx; i--, ii += 2)
					for (j = rry - 1, jj = dy; j >= ry; j--, jj += 2)
					{
						pixel = (double)(imag[ii][jj] + imag[ii + 1][jj] + imag[ii][jj + 1] + imag[ii + 1][jj + 1]) / 4.0;
						imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
					}
				break;
			case R_VERTICAL:
				for (i = rx, ii = dx; i < rrx; i++, ii += 2)
					for (j = rry - 1, jj = dy; j >= ry; j--, jj += 2)
					{
						pixel = (double)(imag[ii][jj] + imag[ii + 1][jj] + imag[ii][jj + 1] + imag[ii + 1][jj + 1]) / 4.0;
						imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
					}
				break;
			case R_HORIZONTAL:
				for (i = rrx - 1, ii = dx; i >= rx; i--, ii += 2)
					for (j = ry, jj = dy; j < rry; j++, jj += 2)
					{
						pixel = (double)(imag[ii][jj] + imag[ii + 1][jj] + imag[ii][jj + 1] + imag[ii + 1][jj + 1]) / 4.0;
						imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
					}
				break;
			case F_DIAGONAL:
				for (j = ry, ii = dx; j < rry; j++, ii += 2)
					for (i = rx, jj = dy; i < rrx; i++, jj += 2)
					{
						pixel = (double)(imag[ii][jj] + imag[ii + 1][jj] + imag[ii][jj + 1] + imag[ii + 1][jj + 1]) / 4.0;
						imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
					}
				break;
			case S_DIAGONAL:
				for (j = rry - 1, ii = dx; j >= ry; j--, ii += 2)
					for (i = rrx - 1, jj = dy; i >= rx; i--, jj += 2)
					{
						pixel = (double)(imag[ii][jj] + imag[ii + 1][jj] + imag[ii][jj + 1] + imag[ii + 1][jj + 1]) / 4.0;
						imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
					}
				break;
			}
		}

		swap(image1, image, PIXEL **);
	}
	printf("\n");
}

void piramidal_decoding(int level)
{
	int rx, ry, rrx, rry, dx, dy;
	int i, j, ii, jj;
	PIXEL **imag, **imag1;
	double pixel;

	if (level < 1)
		return;

	zooming(2.0); /* Increase resolution */

	imag = image;
	imag1 = image1;

	printf(" Increasing resolution... \r");
	fflush(stdout);
	trans = &fractal_code;

	while (trans->next != NULL)
	{
		trans = trans->next;

		rx = (int)rint(trans->rx);
		ry = (int)rint(trans->ry);
		dx = (int)rint(trans->dx);
		dy = (int)rint(trans->dy);
		rrx = (int)rint(trans->rrx);
		rry = (int)rint(trans->rry);

		switch (trans->sym_op)
		{
		case IDENTITY:
			for (i = rx, ii = dx >> 1; i < rrx; i++, ii++)
				for (j = ry, jj = dy >> 1; j < rry; j++, jj++)
				{
					pixel = (double)imag[ii][jj];
					imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
				}
			break;
		case R_ROTATE90:
			for (j = rry - 1, ii = dx >> 1; j >= ry; j--, ii++)
				for (i = rx, jj = dy >> 1; i < rrx; i++, jj++)
				{
					pixel = (double)imag[ii][jj];
					imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
				}
			break;
		case L_ROTATE90:
			for (j = ry, ii = dx >> 1; j < rry; j++, ii++)
				for (i = rrx - 1, jj = dy >> 1; i >= rx; i--, jj++)
				{
					pixel = (double)imag[ii][jj];
					imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
				}
			break;
		case ROTATE180:
			for (i = rrx - 1, ii = dx >> 1; i >= rx; i--, ii++)
				for (j = rry - 1, jj = dy >> 1; j >= ry; j--, jj++)
				{
					pixel = (double)imag[ii][jj];
					imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
				}
			break;
		case R_VERTICAL:
			for (i = rx, ii = dx >> 1; i < rrx; i++, ii++)
				for (j = rry - 1, jj = dy >> 1; j >= ry; j--, jj++)
				{
					pixel = (double)imag[ii][jj];
					imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
				}
			break;
		case R_HORIZONTAL:
			for (i = rrx - 1, ii = dx >> 1; i >= rx; i--, ii++)
				for (j = ry, jj = dy >> 1; j < rry; j++, jj++)
				{
					pixel = (double)imag[ii][jj];
					imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
				}
			break;
		case F_DIAGONAL:
			for (j = ry, ii = dx >> 1; j < rry; j++, ii++)
				for (i = rx, jj = dy >> 1; i < rrx; i++, jj++)
				{
					pixel = (double)imag[ii][jj];
					imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
				}
			break;
		case S_DIAGONAL:
			for (j = rry - 1, ii = dx >> 1; j >= ry; j--, ii++)
				for (i = rrx - 1, jj = dy >> 1; i >= rx; i--, jj++)
				{
					pixel = (double)imag[ii][jj];
					imag1[i][j] = bound(0.5 + pixel * trans->alfa + trans->beta);
				}
			break;
		}
	}
	swap(image1, image, PIXEL **);

	if (level > 1)
	{
		if (quality)
			iterative_decoding(0, 2, 1.0);
		piramidal_decoding(level - 1);
	}
}

void smooth_image()
{
	double pixel1, pixel2;
	int i, j;
	int w1, w2;
	int rx, ry, rrx, rry;

	printf("\n Postprocessing ...");
	fflush(stdout);

	trans = &fractal_code;
	while (trans->next != NULL)
	{
		trans = trans->next;

		rx = (int)rint(trans->rx);
		ry = (int)rint(trans->ry);
		rrx = (int)rint(trans->rrx);
		rry = (int)rint(trans->rry);

		if (rx == 0 || ry == 0 || (int)trans->size == 1)
			continue;

		if (trans->size == min_size)
		{
			w1 = 5;
			w2 = 1;
		}
		else
		{
			w1 = 2;
			w2 = 1;
		}

		for (i = rx; i < rrx; ++i)
		{
			pixel1 = image[i][ry];
			pixel2 = image[i][ry - 1];
			image[i][ry] = (w1 * pixel1 + w2 * pixel2) / (w1 + w2);
			image[i][ry - 1] = (w2 * pixel1 + w1 * pixel2) / (w1 + w2);
		}

		for (j = ry; j < rry; ++j)
		{
			pixel1 = image[rx][j];
			pixel2 = image[rx - 1][j];
			image[rx][j] = (w1 * pixel1 + w2 * pixel2) / (w1 + w2);
			image[rx - 1][j] = (w2 * pixel1 + w1 * pixel2) / (w1 + w2);
		}
	}
	printf("done \n");
	fflush(stdout);
}

double calc_PSNR(char *img1_path, char *img2_path)
{
    FILE *in1, *in2;

    if ((in1 = fopen(img1_path, "rb")) == NULL || (in2 = fopen(img2_path, "rb")) == NULL)
        fatal("Can't open input file !\n");

    PIXEL **img1, **img2;

    matrix_allocate(img1, image_width, image_height, PIXEL);
    matrix_allocate(img2, image_width, image_height, PIXEL);

    int i, j;

    i = fread(img1[0], sizeof(PIXEL), image_width * image_height, in1);
    j = fread(img2[0], sizeof(PIXEL), image_width * image_height, in2);

    if (i < image_width * image_height || j < image_width * image_height)
        fatal("Not enough input data in the input file !\n");

    fclose(in1);
    fclose(in2);

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

    int i, j;

    i = fread(img1[0], sizeof(PIXEL), image_width * image_height, in1);
    j = fread(img2[0], sizeof(PIXEL), image_width * image_height, in2);

    if (i < image_width * image_height || j < image_width * image_height)
        fatal("Not enough input data in the input file !\n");

    fclose(in1);
    fclose(in2);

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
