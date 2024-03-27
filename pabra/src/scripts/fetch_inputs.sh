#!/bin/bash

srcdir=~/GGG/input-FFS-CN50
tgtdir=./input-FFS-CN50

plist10="10.4.10.11 10.4.10.12 10.4.10.13 10.4.10.14 10.4.10.15 10.4.10.16 10.4.10.17 10.4.10.18 10.4.10.19 10.4.10.20 "
plist9=" 9.4.10.11  9.4.10.12  9.4.10.13  9.4.10.14  9.4.10.15  9.4.10.16  9.4.10.17  9.4.10.18  9.4.10.19  9.4.10.20  "


for p in $plist9 $plist10; do

	idir=$tgtdir/input.$p
	mkdirhier $idir

	echo 0	>> $idir/temp
	echo 23 >> $idir/temp

	cp -v $srcdir/*$p* $idir

	mv -v $idir/inputAI.$p.end.mod $idir/inputAI.end

	echo Writing $idir/runfile.$p
	replace _PARS_ $p < fetch_inputs_runfile.template > $idir/runfile.$p
done
