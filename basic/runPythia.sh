#!/bin/bash
# 1st argument: number of events to generate
# 2nd argument: settings file to use
if [ -z "$SLURM_ARRAY_TASK_ID" ] 
  then I=0
  else
  I=$SLURM_ARRAY_TASK_ID
fi
nEvents=$1
settings=$2
mkdir $I
cd $I
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST  > ../log$I.out
cp ../$settings .
../oniaVsMult $nEvents $settings out $I
mv ../log$I.* .
