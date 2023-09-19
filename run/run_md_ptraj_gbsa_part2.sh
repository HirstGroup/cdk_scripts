#!/bin/bash

set -e

# default values
repeat=NO
test=NO

while [[ $# > 0 ]]
do
key="$1"

case $key in
    -c|--complex)
    complex="$2"
    shift
    ;;
    -r|--repeat)
    repeat="$2"
    shift
    ;;
    --test)
    test="$2"
    shift
    ;;
    *)
    echo "Unknown argument: $1"
    exit 1
    ;;
esac
shift
done

#echo "complex = " $complex "repeat = " $repeat

jobid=$(sbatch --parsable ~/scripts/1gpu.sh ~/cdk_scripts/standard_md2.sh -c $complex -r "$repeat" --test $test)

echo "Submitted batch job $jobid"

sbatch --dependency=afterok:$jobid ~/scripts/1cpu.sh ~/cdk_scripts/run/run_ptraj_gbsa.sh -c $complex -p 2 -r "$repeat"