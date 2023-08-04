complex=$1
repeat=${2:-""} # optional parameter to run repeats, e.g. _1 _2 etc

python ~/cdk_scripts/ptraj.py -i $complex -r "$repeat" --delete
python ~/cdk_scripts/gbsa.py -i $complex -r "$repeat"