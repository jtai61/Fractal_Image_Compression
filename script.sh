#!/bin/bash

# clean bin folder
rm -f bin/enc
rm -f bin/dec
rm -f bin/obj/*.o

# clean ifs folder
rm -f ifs/*.ifs

# clean img/decompressed folder
rm -f img/decompressed/*.pgm
rm -f img/decompressed/*.raw

# # encode image to ifs
# time ./bin/enc -F
# time ./bin/enc -X
# time ./bin/enc -C
# time ./bin/enc -S
# time ./bin/enc -Z
# time ./bin/enc -Y

# # decode ifs to image
# time ./bin/dec -i -r -n 3

