#!/bin/bash

######################################################################
# MASTER START SCRIPT FOR RUNNING pabra WRAPPED AROUND GGG simulator #
######################################################################
#
# by T.R. Sokolowski 2012-2020
#

##############
# ATTENTION! # The initialization of the tree must be adapted to the system size in the startscript!
##############

##### Std. parameter flags for GGG code (4th integer varies weak repression unbinding rate)
##### Flag 10.5.10.10 : weak repression = strong repression; flag 10.5.10.20 : no weak repression at all
parlist="10.5.10.10 10.5.10.11 10.5.10.12 10.5.10.13 10.5.10.14 10.5.10.15 10.5.10.16 10.5.10.17 10.5.10.18 10.5.10.19 10.5.10.20"
parlist="10.5.10.13"

##### Time parameters for pabra
ntraj=10 			# no. of independent starts at the root of the branching tree
simtime=$((24000*$ntraj)) 	# total simulated time [s]
bintime=300 			# size of time bins [s]

##### Cluster parameters
#nodes="NODENAME"
#nopt="-l nodes=$nodes"
#nopt=""

##### Set label string and directories
label="CN15_unp_sym_mw1e-4_pf0.25_${simtime}s_${ntraj}traj_SIC_eqS"
jobsetname=FFS_${label}

date="$(date +"%m-%d_%H-%M-%S" )"
datedir="Pabra+GGG_$date"
echo "datedir = $datedir ..."

comment="unp"

##### Master directory / path
masterdir="/home/sokolowski/Pabra+GGG" # Enter where the code is located on your machine here
                                       # Make sure to use a full path name

##### Choose directory for equilibrated initial conditions / input files
#inputbase=${masterdir}/input-FFS-CN15-SIC-equalStripes-pinned/
inputbase=${masterdir}/input-FFS-CN15-SIC-equalStripes-unpinned/

##### Set other important directories
outputbase=${masterdir}/RUN/
targetbase=${masterdir}/DONE/   # the directory to which everything is moved finally
GGG_src_dir=${masterdir}/GGG/

##### Choose the pabra start scripts (defines which reaction coordinates to use)
#scripts="tpropext_tomek_1coord.py tpropext_tomek_2coords.py tpropext_tomek_2coords_orth.py"
#scripts="tpropext_tomek_1coord.py"
scripts="tpropext_tomek_1coord_sym.py"
#scripts="tpropext_tomek_1coord_prom.py"

##### Compile source code with cluster option of the makefile
echo "Compiling ..."
pushd $GGG_src_dir
make cluster
popd

##### Start all chosen scripts with all parameters chosen
for script in $scripts; do
for pars in $parlist; do

	echo "Defining directories ..."
	randomint=$(echo $RANDOM)
	randomdir="${datedir}_sim${randomint}"
	jobdir="${masterdir}/pabra/src/testruns/"
	inputdir="${inputbase}/input.$pars/"
	outputdir="${outputbase}/${randomdir}/"
	targetdir="${targetbase}/${datedir}/sim_${comment}_${pars}/"
	echo "  jobdir    = $jobdir"
	echo "  inputdir  = $inputdir"
	echo "  outputdir = $outputdir"
	echo "  targetdir = $targetdir"

	jobcmd="./start27.sh $script $inputdir $outputdir $pars $simtime $ntraj $bintime $comment"
	jobname="${jobsetname}_${pars}"

	echo "Preparing input data directories ..."
	cp -v $GGG_src_dir/gillespie.X64 ${inputdir}/gillespie.X

	echo "Creating target directory ..."
	mkdirhier ${targetdir}

	echo "Constructing script for job $jobname ..."
	echo "jobcmd = $jobcmd"
	echo

	echo "Creating output directory ${outputdir}"
	mkdirhier ${outputdir}

	echo "Copying input files"
	cp -v ${inputdir}/* ${outputdir}

	# Enter job directory
	pushd ${jobdir}

	# Run job command with parameters and redirect output
	export $PYTHONPATH
	${jobcmd}

	# Return to previous directory
	popd

	# Copy final output to targetdir
	cp -rv ${outputdir}/* ${targetdir}
done
done


echo "*** DONE! ***"

