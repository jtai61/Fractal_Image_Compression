#!/bin/bash

function clean()
{
	# clean bin folder
	rm -rf bin/

	# clean ifs folder
	rm -f ifs/*.ifs

	# clean img/decompressed folder
	rm -f img/decompressed/*.pgm
	rm -f img/decompressed/*.raw
}


clean

exit 0
