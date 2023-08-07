#!/bin/bash

set -e

# default values
part=NO
repeat=NO
sbatch=YES
test=NO

while [[ $# > 0 ]]
do
key="$1"

case $key in
    -c|--complex)
    complex="$2"
    shift
    ;;
    -p|--part)
    part="$2"
    shift
    ;;
    -r|--repeat)
    repeat="$2"
    shift
    ;;
    --sbatch)
    sbatch="$2"
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

echo "complex = " $complex "part = " $part "repeat = " $repeat "test = " $test

sbatch_cmd="sbatch --parsable ~/scripts/1gpu.sh"
if [ $sbatch = "NO" ]; then
	sbatch_cmd=""
fi

if [ $part = "NO" ]; then
	jobid=$(${sbatch_cmd}bash ~/cdk_scripts/standard_md.sh -c $complex -r "$repeat" --test $test)
else
	jobid=$(${sbatch_cmd}bash ~/cdk_scripts/standard_md${part}.sh -c $complex -r "$repeat" --test $test)
fi

echo "Running batch job $jobid"

sbatch_cmd="sbatch --dependency=afterok:$jobid ~/scripts/1cpu.sh bash "
if [ $sbatch = "NO" ]; then
	sbatch_cmd=""
fi

${sbatch_cmd}bash ~/cdk_scripts/run/run_ptraj_gbsa.sh -c $complex -p "$part" -r "$repeat"