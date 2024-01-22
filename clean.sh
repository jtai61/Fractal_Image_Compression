#!/bin/bash

function clean()
{
	# delete bin and obj folder
	rm -rf bin/
	rm -rf obj/

	# clean ifs folder
	rm -f ifs/*.ifs

	# clean decompressed folder
	rm -f img/decompressed/*.pgm
	rm -f img/decompressed/*.raw
}


clean

exit 0
