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
    -t|--time)
    time="$2"
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

if [$repeat -eq ""]; then
	python ~/cdk_scripts/ptraj.py -i $complex -t "equi$part" --delete
	python ~/cdk_scripts/gbsa.py -i $complex -p "$part"
else
	python ~/cdk_scripts/ptraj.py -i $complex -t "equi$part" --delete
	python ~/cdk_scripts/gbsa.py -i $complex -p "$part"
fi