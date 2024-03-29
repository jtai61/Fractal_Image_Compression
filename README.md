# DESCRIPTION

This program is aimed to compare several speed-up 
techniques in fractal image coding. In particular six methods
have been implemented, chosen from both  classification  and 
feature vectors approaches. The methods are: Fisher, Hurtgen, 
MassCenter, Saupe, Saupe-Fisher and MassCenter-Saupe. Please
refer to the MANUAL for further details. 
I tried to keep it simple, many solutions have been adopted
for their simplicity and not because they are the best. 
Some fragments of code are repeted several times. This is because
I tried to make the various speed-up methods indipendent each
other in such a way you can see each method as a plug-in. So
if you want to try a new method you need just to write the code
that implement it and to insert it in the coder. In particular
you need to write only two functions.

Other features supported:
- Entropy based split decision function
- Variance based split decision function
- Adaptive splitting thresholds
- Output quadtree partition in pgm format
- Iterative decoding
- Piramidal decoding
- Fractal zooming

# ENCODER

Encode a grayscale image using a quadtree based fractal image coder.
The input image is either in raw or in pgm format.

## USAGE
        enc [ -options ] [ inputfile [ outputfile ] ]

## OPTIONS

### -F
- The Fisher speed-up method is used [1].

### -X
- A derivative of the Hurtgen [2] speed-up method is used. It 
works as follow. A block is partitionned in its four quadrants 
then the mean of each quadrant is computed. To each quadrant 
is assigned a bit which is 1 if its mean is above the overall
mean and 0 otherwise. This way a block is associated with one
of 16 possible classes. Indeed the class 1111 is always empty
so there are only 15 classes. To further characterize the 
blocks, for each of the above classes, 24 subclasses are 
considered by using the variance method used in the Fisher
classification scheme yelding 360 classes in all.

### -C
- The Mass center speed-up method [3] is used.

### -S
- The Saupe speed-up method [4] is used.

### -Z
- The Saupe speed-up method is used along with the isometries 
computation technique of Fisher.

### -Y
- The Saupe speed-up method is used combined with the Mass center
method. 
- NOTE: with this method small values of the -n flag should
be used. 

### -r
- Rms threshold used for driving the partition process by the
standard rms based split decision function.

### -e
- Use an entropy based split decision function for partitionning 
the image with the indicated threshold [5].
- NOTE: in order to 
use entropy for driving the partitionning process, the rms split
decision function should be disabled. This can be achieved by
giving to the flag -r a big value.

### -v
- Use a variance based split decision function for partitionning 
the image with the indicated threshold [5].
- NOTE: in order to 
use variance for driving the partitionning process, the rms split
decision function should be disabled. This can be achieved by
giving to the flag -r a big value.

### -a
- Adaptive threshold factor used by splitting decision 
functions [5]. Both rms and variance functions use the 
following adaptive threshlold relation: 
T(i) = adaptive_factor * T(i-1) where T(i) is the threshold on 
the level i of the quadtree. For entropy instead the following 
threshold relation is used: T(i) = T(i-1) + adaptive_factor / i.
It seems that the best adaptive values are 2.0 and 4.0 for 
respectively rms and variance functions. For the entropy function
the best value seems to be in the range [sqrt(2.0),2.0] depending 
on the number of levels of the quadtree.
It is interesting to see how adaptivity yelds better quality
decoded images by attenuating the effect of blockiness.

### -d
- Domain step used for building the codebook. Only even values
are supported.

### -W
- Image width

### -H
- Image height

### -m
- Min range size in the quadtree.

### -M
- Max range size in the quadtree.

### -l
- Number of neighbours requested in a kd-tree query

### -p
- Tolerance used in a kd-tree query

### -f
- Full 1^ class search used in both the Fisher and Hurtgen method.

### -s
- Full 2^ class search used in both the Fisher and Hurtgen method.

### -k
- Shrunk factor in MC method. If greater than 1 a block is shrunken
before computing the mass center.

### -c
- This flag is related to the number of Saupe features. It has two
meanings. If it is less than 4 then the block will be shrunken by 
a factor of 2^flag_value and then the features will be extracted 
from the resulting block. Otherwise the flag indicates directly the
number of features to be used, the shrunk factor will be computed
automatically, only power of 4 value are accepted.

### -n
- Number of classes used in the Mass center method.

### -A
- Number of bits used for quantizing the scaling factor of a fractal
transformation.

### -B
- Number of bits used for quantizing the offset of a fractal
transformation.

### -y
- Max allowed value for the scaling factor. NOTE: only positive
values are supported for the scaling factor. This avoid to search
twice the codebook.

### -Q
- Output the quadtree partition as image in pgm format.

### -z
- Force the current range block to be considered flat. And hence
approssimated with the scaling factor equal to 0.

### -h
- Help 


# DECODER

Decode a fractal code using either an iterative algorithm or 
a piramidal [6] algorithm.

## USAGE
        dec [ -options ] [ inputfile [ outputfile ] ]

## OPTIONS

### -n
- Number of iterations used to decode the image. 

### -p
- Performs a postprocessing on the decoded image.

### -q
- Better quality decoding when piramidal decoding is used. Combines 
piramidal with iterative decoding. Advised when zooming the image.

### -i
- Iterative decoding   

### -r
- Output image in raw format

### -d
- Display the decoded image by piping it on xv

### -z
- Zoom factor. Decode the image at a resolution different from
the original one.

### -h
- Help



# REFERENCES

[1] Y. Fisher, ``Fractal Image Compression--Theory and Application'',
        Springer-Verlag, New York, 1994.

[2] B. Hurtgen, C. Stiller, ``Fast Hierarchical Codebook Search For
        Fractal Coding of Still Image'', EOS/SPIE Visual Communication
	and PACS for medical applications '93, Berlin 1993.

[3] M. Polvere, M. Nappi, ``A Feature Vector Tecnhique For Fast Fractal
        Image Coding'', technical report University of Salerno 1998.

[4] D. Saupe, ``Fractal image compression by multi-dimensional nearest 
        neighbor search'', Proceedings DCC'95 Data Compression Conference, 
	J. A. Storer and M. Cohn (eds.), IEEE Comp. Soc. Press, March 1995.

[5] R. Distasi, M. Polvere, M. Nappi, ``Split Decision Functions in 
        Fractal Image Coding'', Electronics Letters 34,8, vol. 34, no. 8, 
	pp. 751--753, April 1998.

[6] Z. Baharav, D. Malah, E. Karnini, ``Hierarchical Interpretation
        of Fractal Image Coding and its Application to Fast Decoding'',
	In Inlt. Conference on Digital Signal Processing, Cyprus July 1993.

