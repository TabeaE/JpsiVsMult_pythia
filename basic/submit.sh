#!/bin/bash
# example use:
# ./submit.sh outputFolder 1000 1000000 settings.cmnd 
# first argument: output folder
# 2nd argument: number of jobs (size of job array)
# 3rd argument: number of events per job
# 4th argument: settings file to use
outDir=$1
nJobs=$2
nEvents=$3
settings=$4
mkdir trees/$outDir
cd trees/$outDir
cp ../../oniaVsMult .
cp ../../oniaVsMult.cc .
cp ../../runPythia.sh .
cp ../../$settings .
sbatch --mem-per-cpu=3072 --array=1-$nJobs -p long -o log%a.out -e log%a.err runPythia.sh $nEvents $settings
