complex=$1
repeat=${2:-""} # optional parameter to run repeats, e.g. _1 _2 etc

echo "complex = " $complex "repeat = " $repeat

jobid=$(sbatch --parsable ~/scripts/1gpu.sh ~/cdk_scripts/standard_md.sh $complex "$repeat")

jobid=$(sbatch --dependency=afterok:$jobid --parsable ~/scripts/1cpu.sh ~/cdk_scripts/ptraj.py -i $complex -r "$repeat")

jobid=$(sbatch --dependency=afterok:$jobid --parsable ~/scripts/1cpu.sh ~/cdk_scripts/gbsa.py -i $complex -r "$repeat")





