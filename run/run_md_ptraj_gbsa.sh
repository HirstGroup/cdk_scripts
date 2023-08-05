#!/bin/bash

complex=$1
repeat=${2:-""} # optional parameter to run repeats, e.g. _1 _2 etc

echo "complex = " $complex "repeat = " $repeat

jobid=$(sbatch --parsable ~/scripts/1gpu.sh bash ~/cdk_scripts/standard_md.sh $complex "$repeat")

sbatch --dependency=afterok:$jobid ~/scripts/1cpu.sh bash ~/cdk_scripts/run/run_ptraj_gbsa.sh $complex "$repeat"





