#!/bin/bash

dumpscript='dump_tree.py'
treename='tpropext.h5'
histtime=15000

#task='segments'
#task='paths'
task='histograms'
#task='correlations'

pythonpath='/home/thomas/FFS/pabra/src'

ziplist="gillespie* inputAI* reactions* runfile* observables* temp*"

for d in $*; do

	cp -v $dumpscript $d

	pushd $d

		if [[ "$task" == "paths"  ]]; then
			rm -r ./paths
		fi

		if [[ "$task" == "segments"  ]]; then
			rm -r ./segments
		fi

		# FOR ALL TASKS DO:
		PYTHONPATH=$pythonpath python2.7 $dumpscript $treename $task $histtime

		if [[ "$task" == "paths"  ]]; then
			gnuplot ./paths/plot_paths.gnu
		fi

		if [[ "$task" == "segments"  ]]; then
			gnuplot ./segments/plot_segments.gnu
			gnuplot ./segments/plot_phasespace.gnu
		fi

		if [[ "$task" == "histograms"  ]]; then
			gnuplot ./trunkdata/plot_hist_rc.gnu
			gnuplot ./trunkdata/plot_hist_ps.gnu
			gnuplot ./trunkdata/plot_hist_rc_cnt.gnu
			gnuplot ./trunkdata/plot_hist_ps_cnt.gnu
			gnuplot ./trunkdata/plot_running_avg.gnu
			gnuplot ./trunkdata/plot_basin_occupancies.gnu
		fi

		if [[ "$task" == "correlations"  ]]; then
			gnuplot ./corr_func/plot_corr_func.gnu
		fi

#		echo Zipping and removing schools ...
#		zip -r schools.zip s*x0
#		rm -r s*x0

		echo Zipping input files ...
		for z in $ziplist; do

			if [ -f $z ]; then
				zip -r input.zip $z
				rm $z
			fi
		done
	popd

done
