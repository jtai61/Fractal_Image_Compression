#!/bin/bash

function clean()
{
	# delete bin and obj folder
	rm -rf build/

	# clean ifs folder
	rm -f ifs/*.ifs

	# clean decompressed folder
	rm -f img/decompressed/*.pgm
	rm -f img/decompressed/*.raw
}

function compile()
{
	# create build folder if not exist
	if [ ! -d "build/" ]; then
		mkdir build
		mkdir build/obj
	fi

	# compile
	make
}

function encode()
{
	if [ -f "./build/enc" ]; then
		# encode raw to ifs
		./build/enc $1
	else
		echo "enc file does not exist !"
		exit 1
	fi
}

function decode()
{
	if [ -f "./build/dec" ]; then
		# decode ifs to raw
		./build/dec -i -r
	else
		echo "dec file does not exist !"
		exit 1
	fi
}

# parse args
while [[ $# -ge 1 ]]; do
	mode="$1"
	case $mode in
		-c )
			clean
			;;
		-m )
			clean
			compile
			;;
		-e )
			encode "$2"
			shift
			;;
		-d )
			decode
			;;
	esac
	shift
done
