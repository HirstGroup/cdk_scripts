#!/bin/bash

complex=$1
repeat=${2:-""} # optional parameter to run repeats, e.g. _1 _2 etc
part=${3:-""} # optional parameter to run second part of equi with standard_md2

python ~/cdk_scripts/ptraj.py -i $complex -r "$repeat" -t "equi$part" --delete
python ~/cdk_scripts/gbsa.py -i $complex -p "$part" -r "$repeat"