#!/bin/bash

dumpscript='dump_tree.py'
treename='tpropext.h5'
pythonpath='/home/sokolowski-shared/data/TEST_FFS/pabra/src'

# Specify simulation time interval and bins.
# Make sure that histtime/nbins divides the simulation time bin size
# for correct working of the data binning.
histtime=24000 # 6.66h data
nTbins=40

nHbins=100	# phase space histogram bins
tLog=300	# the logging interval of the FFS sim. in seconds
tUpdate=60	# the time between two FFS restarts

plot_only=0	# whether to regenerate the plots only

# Specify which analysis routine to run.
# Can be any combination from keywords: segments paths histograms correlations
tasklist='histograms'

ziplist="gillespie* inputAI* reactions* runfile* observables* temp*"

for d in $*; do
for task in $tasklist; do

	cp -v $dumpscript $d

	pushd $d

		if [[ "$task" == "paths"  ]]; then
			rm -r ./paths
		fi

		if [[ "$task" == "segments"  ]]; then
			rm -r ./segments
		fi

		# FOR ALL TASKS DO:
		PYTHONPATH=$pythonpath python2.7 $dumpscript $treename $task $histtime $nTbins $nHbins $tLog $tUpdate $plot_only

		if [[ "$task" == "paths"  ]]; then
			gnuplot ./paths/plot_paths.gnu
			mv -v paths.ps paths.pdf paths/
			cp ../plot_snapshots.sh ./paths/
			pushd paths

				#./plot_snapshots.sh
				for psf in *.ps; do ps2pdf $psf; done
			popd
		fi

		if [[ "$task" == "segments"  ]]; then
			gnuplot ./segments/plot_segments.gnu
			gnuplot ./segments/plot_phasespace.gnu
			mv -v segments.ps phasespace.ps segments/
			pushd segments

				for psf in *.ps; do ps2pdf $psf; done
			popd
		fi

		if [[ "$task" == "histograms"  ]]; then

			for pf in plot_hist_rc plot_hist_ps plot_hist_rc_cnt plot_hist_ps_cnt plot_running_avg plot_basin_occupancies; do

				gnuplot ./trunkdata_$nHbins/${pf}.gnu
			done
			for psf in *.ps; do ps2pdf $psf; done
			mv -v *.ps *.pdf trunkdata_$nHbins
		fi

		if [[ "$task" == "correlations"  ]]; then
			gnuplot ./corr_func/plot_corr_func.gnu
		fi

		echo "Zipping input files for archiving ..."
		for z in $ziplist; do

			if [ -f $z ]; then
				zip -r input.zip $z
				rm $z
			fi
		done
	popd


done # tasks
done # dirs

echo "Done!"
