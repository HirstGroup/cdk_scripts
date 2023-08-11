#!/bin/bash

set -e

# default values
ignore_errors=NO
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
    --ignore_errors)
    ignore_errors="$2"
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

echo $ignore

python ~/cdk_scripts/ptraj.py -i $complex -r "$repeat" -p $part --delete --ignore_errors ${ignore_errors} --test $test
python ~/cdk_scripts/gbsa.py -i $complex -p "$part" -r "$repeat" --test $test