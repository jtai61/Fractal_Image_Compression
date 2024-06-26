#include "base.h"

void fatal(const char *s)
{
	printf("%s\n", s);
	exit(-1);
}

void help_enc()
{
	printf(
		"\n Usage: enc [-options] [inputfile [outputfile]]"
		"\n Supported options:\n"
		" --------------------------------------------------------------------------\n"
		" -F   Fisher method          (Off)   -X   Hurtgen method              (Off)\n"
		" -Z   Saupe-Fisher method    (Off)   -S   Saupe method                (Off)\n"
		" -Y   Saupe-MC method        (Off)   -C   MassCenter method           (On) \n"
		" -T   Tai method             (Off)	  -D   Baseline method			   (Off)\n"
		" -e # Entropy threshold      (8.0)   -d # Domain step                 (4)  \n"
		" -r # Rms threshold          (8.0)   -a # Adaptive factor             (1.0)\n"
		" -v # Variance threshold     (Inf)   -W # Image width                 (512)\n"
		" -m # Min range size         (4)     -H # Image height                (512)\n"
		" -M # Max range size         (16)    -l # Matches in a kd-tree query  (50) \n"
		" -f   Full 1^ class search   (Off)   -p # Epsilon in a kd-tree query  (2.0)\n"
		" -s   Full 2^ class search   (Off)   -k # Shrunk factor in MC method  (1.0)\n"
		" -c # Num of Saupe features  (16)    -n # Number of MC classes        (50) \n"
		" -A # Num of bits for alfa   (4)     -B # Number of bits for beta     (7)  \n"
		" -y # Max value of alfa      (1.0)   -Q   Output a quadtree image     (Off)\n"
		" -z # Force alfa to be 0     (0)     -h   Display this help \n\n"
		" Supported image format: raw and pgm\n"
		" Default input file : img/uncompressed/Lena_256.raw\n"
		" Default output file: ifs/Lena_256_enc.ifs\n\n");
}

void help_dec()
{
	printf("\n Usage: dec [-options] [inputfile [outputfile]]"
		   "\n Supported options:\n"
		   " --------------------------------------\n"
		   " -n # Num of iterations           (10) \n"
		   " -p   Postprocessing image        (Off)\n"
		   " -q   Better quality decode       (Off)\n"
		   " -i   Iterative decoding          (Off)\n"
		   " -r   Output image in raw format  (Off)\n"
		   " -d   Display image with xv       (Off)\n"
		   " -z # Zoom factor                 (1.0)\n"
		   " -h   Display this help \n\n"
		   " Default input file : ifs/Lena_256_enc.ifs\n"
		   " Default output file: img/decompressed/Lena_256_dec.pgm\n\n");
}

void getopt_enc(int argc, char **argv)
{
	int i, num;

	filein[0] = 1;
	fileout[0] = 1;

	for (i = 1; i < argc; ++i)
	{
		if (argv[i][0] != '-')
		{
			if (filein[0] == 1)
				strcpy(filein, argv[i]);
			else if (fileout[0] == 1)
				strcpy(fileout, argv[i]);
			else
				;
		}
		else
		{
			if (strlen(argv[i]) == 1)
				break;
			switch (argv[i][1])
			{
			case 'D':
				method = Baseline;
				break;
			case 'F':
				method = Fisher;
				break;
			case 'X':
				method = Hurtgen;
				break;
			case 'Z':
				method = SaupeFisher;
				break;
			case 'S':
				method = Saupe;
				break;
			case 'C':
				method = MassCenter;
				break;
			case 'Y':
				method = McSaupe;
				break;
			case 'T':
				method = Tai;
				break;
			case 'r':
				T_RMS = atof(argv[++i]);
				break;
			case 'e':
				T_ENT = atof(argv[++i]);
				break;
			case 'v':
				T_VAR = atof(argv[++i]);
				break;
			case 'a':
				adapt = atof(argv[++i]);
				if (adapt <= 0.0)
					fatal("\n -a flag not greater than 0.0");
				break;
			case 'H':
				image_height = atoi(argv[++i]);
				break;
			case 'W':
				image_width = atoi(argv[++i]);
				break;
			case 'd':
				SHIFT = atoi(argv[++i]);
				if (SHIFT < 2 || SHIFT > 32 || SHIFT % 2)
					fatal("\n -d flag not in the range [2,32] or not even");
				break;
			case 'm':
				min_size = atoi(argv[++i]);
				if (min_size < 2 || min_size > 64)
					fatal("\n -m flag not in the range [2,64]");
				num = min_size;
				while (num > 2)
				{
					if ((num % 2))
						fatal("\n -m flag not a power of two");
					num /= 2;
				}
				if (num != 2)
					fatal("\n -m flag not a power of two");
				break;
			case 'M':
				max_size = atoi(argv[++i]);
				if (max_size < 2 || max_size > 64)
					fatal("\n -M flag not in the range [2,64]");
				num = max_size;
				while (num > 2)
				{
					if ((num % 2))
						fatal("\n -m flag not a power of two");
					num /= 2;
				}
				if (num != 2)
					fatal("\n -M flag not a power of two");
				break;
			case 'l':
				matches = atoi(argv[++i]);
				if (matches < 1 || matches > MAX_NEIGHBOURS)
					fatal("\n -l flag not in the range [1,MAX_NEIGHBOURS]");
				break;
			case 'p':
				eps = atof(argv[++i]);
				if (eps <= 0.0 || eps > 20.0)
					fatal("\n -p flag not in the range ]0.0,20.0]");
				break;
			case 'f':
				full_first_class = 1;
				break;
			case 's':
				full_second_class = 1;
				break;
			case 'k':
				shrunk_factor_mc = atoi(argv[++i]);
				if (shrunk_factor_mc < 1)
					fatal("\n -k flag less than 1");
				num = shrunk_factor_mc;
				while (num > 1)
				{
					if ((num % 2))
						fatal("\n -k flag not a power of two");
					num /= 2;
				}
				if (num != 1)
					fatal("\n -k flag not a power of two");
				break;
			case 'c':
				n_features = atoi(argv[++i]); /* MaxFeatures */
				if (n_features == 0)
					fatal("\n -c flag not allowed");
				if (n_features > 0 && n_features < 4)
					shrunk_factor_saupe = 1 << n_features;
				else
				{
					num = n_features;
					while (num > 4)
					{
						if ((num % 4))
							fatal("\n -k flag not a power of two");
						num /= 4;
					}
					if (num != 4)
						fatal("\n -c flag not a power of four");
				}
				break;
			case 'n':
				n_p_class = atoi(argv[++i]);
				if (n_p_class <= 0)
					fatal("\n -n flag not greater than 0");
				break;
			case 'A':
				N_BITALFA = atoi(argv[++i]);
				if (N_BITALFA < 1 || N_BITALFA > 15)
					fatal("\n -A flag not in the range [1,15]");
				break;
			case 'B':
				N_BITBETA = atoi(argv[++i]);
				if (N_BITBETA < 1 || N_BITBETA > 15)
					fatal("\n -B flag not in the range [1,15]");
				break;
			case 'y':
				MAX_ALFA = atof(argv[++i]);
				if (MAX_ALFA < 0.0 || MAX_ALFA > 8.0)
					fatal("\n -Y flag not in the range [0.0,5.0]");
				break;
			case 'z':
				zero_threshold = atoi(argv[++i]);
				if (zero_threshold < 0)
					fatal("\n -z flag negative");

				break;
			case 'Q':
				qtree = 1;
				break;
			case 'h':
			default:
				help_enc();
				exit(-1);
			}
		}
	}

	if (filein[0] == 1)
		strcpy(filein, "/home/joe/Fractal_Image_Compression/img/uncompressed/Lena_512.raw");

	if (fileout[0] == 1)
		strcpy(fileout, "/home/joe/Fractal_Image_Compression/ifs/Lena_512_enc.ifs");

	if (min_size > max_size)
		fatal("\n -m flag value greater tham -M flag value");
}

void getopt_dec(int argc, char **argv)
{
	int i;

	filein[0] = 1;
	fileout[0] = 1;

	for (i = 1; i < argc; ++i)
	{
		if (argv[i][0] != '-')
		{
			if (filein[0] == 1)
				strcpy(filein, argv[i]);
			else if (fileout[0] == 1)
				strcpy(fileout, argv[i]);
			else
				;
		}
		else
		{
			if (strlen(argv[i]) == 1)
				break;
			switch (argv[i][1])
			{
			case 'n':
				iterations = atoi(argv[++i]);
				break;
			case 'p':
				postproc = 1;
				break;
			case 'q':
				quality = 1;
				break;
			case 'i':
				piramidal = 0;
				break;
			case 'r':
				raw_format = 1;
				break;
			case 'd':
				display = 1;
				break;
			case 'z':
				zoom = atof(argv[++i]);
				break;
			case 'h':
			default:
				help_dec();
				exit(-1);
			}
		}
	}

	if (filein[0] == 1)
		strcpy(filein, "/home/joe/Fractal_Image_Compression/ifs/Lena_512_enc.ifs");

	if (fileout[0] == 1)
	{
		if (raw_format)
			strcpy(fileout, "/home/joe/Fractal_Image_Compression/img/decompressed/Lena_512_dec.raw");
		else
			strcpy(fileout, "/home/joe/Fractal_Image_Compression/img/decompressed/Lena_512_dec.pgm");
	}
}
