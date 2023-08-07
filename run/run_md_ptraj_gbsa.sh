#!/bin/bash

set -e

# default values
part=NO
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
    -p|--part)
    part="$2"
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

echo "complex = " $complex "part = " $part "repeat = " $repeat "test = " $test

sbatch="sbatch --parsable ~/scripts/1gpu.sh"
if [ $test = "YES" ]; then
	sbatch=""
fi

if [ $part = "NO" ]; then
	jobid=$(${sbatch}bash ~/cdk_scripts/standard_md.sh -c $complex -r "$repeat" --test $test)
else
	jobid=$(${sbatch}bash ~/cdk_scripts/standard_md${part}.sh -c $complex -r "$repeat" --test $test)
fi

echo "Running batch job $jobid"

sbatch="sbatch --dependency=afterok:$jobid ~/scripts/1cpu.sh bash "
if [ $test = "YES" ]; then
	sbatch=""
fi

${sbatch}bash ~/cdk_scripts/run/run_ptraj_gbsa.sh -c $complex -p "$part" -r "$repeat"