#!/bin/bash

set -e

# default values
part=NO
repeat=NO

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
    *)
    echo "Unknown argument: $1"
    exit 1
    ;;
esac
shift
done

python ~/cdk_scripts/ptraj.py -i $complex -r "$repeat" -p $part --delete
python ~/cdk_scripts/gbsa.py -i $complex -p "$part" -r "$repeat"