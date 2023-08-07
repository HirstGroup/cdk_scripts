#!/bin/bash

set -e

# default values
part=""
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
    test=YES
    ;;
    *)
    echo "Unknown argument: $1"
    exit 1
    ;;
esac
shift
done

echo "complex = " $complex "part = " $part "repeat = " $repeat "test = " $test

if [$repeat -eq ""]; then
	jobid=$(sbatch --parsable ~/scripts/1gpu.sh bash ~/cdk_scripts/standard_md${part}.sh -c $complex --test $test)
else
	jobid=$(sbatch --parsable ~/scripts/1gpu.sh bash ~/cdk_scripts/standard_md${part}.sh -c $complex -r "$repeat" --test $test)
fi

if [$repeat -eq ""]; then
	sbatch --dependency=afterok:$jobid ~/scripts/1cpu.sh bash ~/cdk_scripts/run/run_ptraj_gbsa.sh -c $complex -p "$part"
else
	sbatch --dependency=afterok:$jobid ~/scripts/1cpu.sh bash ~/cdk_scripts/run/run_ptraj_gbsa.sh -c $complex -p "$part" -r "$repeat"
fi





