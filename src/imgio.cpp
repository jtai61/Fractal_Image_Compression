#include "base.hpp"

void writeimage_pgm(const char *nome, PIXEL **imag, int width, int height)
{
	FILE *out;
	int i, j;

	if ((out = fopen(nome, "wb")) == NULL)
		fatal("\n Can't open output file");

	printf("\n Writing %s ... ", nome);
	fflush(stdout);
	fprintf(out, "P%d\n%d %d\n255\n", 5, width, height);
	for (i = 0; i < height; i++)
		for (j = 0; j < width; j++)
			fputc((int)imag[i][j], out);

	fclose(out);
	printf("done\n");
}

void writeimage_raw(const char *nome, PIXEL **imag, int width, int height)
{
	FILE *out;
	int i, j;

	if ((out = fopen(nome, "wb")) == NULL)
		fatal("\n Can't open output file");

	printf("\n Writing %s ... ", nome);
	fflush(stdout);
	for (i = 0; i < height; i++)
		for (j = 0; j < width; j++)
			fputc((int)imag[i][j], out);

	printf("done\n");
	fclose(out);
}

void writeimage_pipe(FILE *out, PIXEL **imag, int width, int height)
{
	int i, j;

	fprintf(out, "P%d\n%d %d\n255\n", 5, width, height);
	for (i = 0; i < height; i++)
		for (j = 0; j < width; j++)
			fputc((int)imag[i][j], out);

	fclose(out);
}

int readint(FILE *inp)
{
	int c, n;

	while (1)
	{
		c = getc(inp);
		if (isdigit(c))
		{
			ungetc(c, inp);
			break;
		}
		if (c == '#')
			while (c >= ' ')
				c = getc(inp);
		if (c == EOF)
			return -1;
	}
	fscanf(inp, "%d", &n);
	return n;
}

void readimage_pgm(char *name, int *width, int *height)
{
	FILE *inp;
	int c, c1, i, h, w;

	if ((inp = fopen(name, "rb")) == NULL)
		fatal("\n Can't open input file");

	c = getc(inp);
	c1 = getc(inp);
	if (c != 'P' || (c1 != '5' && c1 != '6'))
	{
		fclose(inp);
		fatal("\n File not in pgm format");
	}

	if (c1 == '6')
		fatal("\n Color images not supported ");

	w = readint(inp);
	h = readint(inp);

	readint(inp);
	/* eat remaining whitespace */
	while ((c = getc(inp)) >= ' ')
		;

	matrix_allocate(image, w, h, PIXEL);

	printf("\n Reading %s (%dx%d) ... ", name, w, h);
	fflush(stdout);

	i = fread(image[0], sizeof(PIXEL), h * w, inp);
	if (i < h * w)
		fatal("\n Not enough input data in the input file");

	printf("done\n");

	*width = w;
	*height = h;

	fclose(inp);
}

void readimage_raw(char *name)
{
	FILE *in;
	int i;

	if ((in = fopen(name, "rb")) == NULL)
		fatal("\n Can't open input file");

	matrix_allocate(image, image_width, image_height, PIXEL);

	printf("\n Reading %s (%dx%d) ... ", name, image_width, image_height);
	fflush(stdout);

	i = fread(image[0], sizeof(PIXEL), image_width * image_height, in);
	if (i < image_width * image_height)
		fatal("\n Not enough input data in the input file");
	fclose(in);
	printf("done\n");
}

int pack(int size, long value, FILE *foutf)
{
	int i;
	static int ptr = 1, sum = 0, num_of_packed_bytes = 0;

	/* size == -1 means we are at the end, so write out what is left */
	if (size == -1 && ptr != 1)
	{
		fputc(sum << (8 - ptr), foutf);
		++num_of_packed_bytes;
		return (0);
	}

	/* size == -2 means we want to know how many bytes we have written */
	if (size == -2)
		return (num_of_packed_bytes);

	for (i = 0; i < size; ++i, ++ptr, value = value >> 1, sum = sum << 1)
	{
		if (value & 1)
			sum |= 1;

		if (ptr == 8)
		{
			fputc(sum, foutf);
			++num_of_packed_bytes;
			sum = 0;
			ptr = 0;
		}
	}
	return (-1);
}

long unpack(int size, FILE *fin)
{
	int i;
	int value = 0;
	static int ptr = 1; /* how many bits are packed in sum so far */
	static int sum;

	/* size == -2 means we initialize things */
	if (size == -2)
	{
		sum = fgetc(fin);
		sum <<= 1;
		return ((long)0);
	}

	/* size == -1 means we want to peek at the next bit without */
	/* advancing the pointer */
	if (size == -1)
		return ((long)((sum & 256) >> 8));

	for (i = 0; i < size; ++i, ++ptr, sum <<= 1)
	{
		if (sum & 256)
			value |= 1 << i;

		if (ptr == 8)
		{
			sum = getc(fin);
			ptr = 0;
		}
	}
	return ((long)value);
}
