#!/bin/bash


dirlist="$1"
tag="$BASH_ARGV"

# Get most recent analysis script
cp -v pabra/src/scripts/dump_tree.py .

# Go through the directories and do
# the analysis
for d in $dirlist*; do

	cp -v dump* $d
	pushd $d

	./dump_all.sh $tag*

	popd

done

