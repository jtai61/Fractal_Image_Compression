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

function start()
{
	if [ -f "./build/enc" ] && [ -f "./build/dec" ]; then

		# # Baseline method
		# echo "--------------------------------------------------Encode--------------------------------------------------"
		# ./build/enc -D
		# echo "--------------------------------------------------Decode--------------------------------------------------"
		# ./build/dec -i -r

		# Fisher method
		echo "--------------------------------------------------Encode--------------------------------------------------"
		./build/enc -F
		echo "--------------------------------------------------Decode--------------------------------------------------"
		./build/dec -i -r

		# Hurtgen method
		echo "--------------------------------------------------Encode--------------------------------------------------"
		./build/enc -X
		echo "--------------------------------------------------Decode--------------------------------------------------"
		./build/dec -i -r

		# Saupe-Fisher method
		echo "--------------------------------------------------Encode--------------------------------------------------"
		./build/enc -Z
		echo "--------------------------------------------------Decode--------------------------------------------------"
		./build/dec -i -r

		# Saupe method
		echo "--------------------------------------------------Encode--------------------------------------------------"
		./build/enc -S
		echo "--------------------------------------------------Decode--------------------------------------------------"
		./build/dec -i -r

		# Saupe-MC method
		echo "--------------------------------------------------Encode--------------------------------------------------"
		./build/enc -Y
		echo "--------------------------------------------------Decode--------------------------------------------------"
		./build/dec -i -r

		# MassCenter method
		echo "--------------------------------------------------Encode--------------------------------------------------"
		./build/enc -C
		echo "--------------------------------------------------Decode--------------------------------------------------"
		./build/dec -i -r

		# Tai method
		echo "--------------------------------------------------Encode--------------------------------------------------"
		./build/enc -T -l 20
		echo "--------------------------------------------------Decode--------------------------------------------------"
		./build/dec -i -r

	else
		echo "enc/dec file does not exist !"
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
		-s )
			start
			;;
	esac
	shift
done
