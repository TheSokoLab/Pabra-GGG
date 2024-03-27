#!/bin/bash


dirlist=$1
tag=$2

# Go through the directories
for d in $dirlist*; do

	cp -v clean_schools $d
	pushd $d

	# Clean in all subdirectories
	for sd in $tag*; do

		echo Cleaning $sd ...
		./clean_schools $sd
	done

	popd

done

