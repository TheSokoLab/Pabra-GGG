#!/bin/bash


dirlist=$(ls -R | grep '/' | cut -d':' -f1)
infile="inputAI*"

for d in $dirlist; do

	echo -n $d : 
	./calc_lambda.sh $d/$infile

done
