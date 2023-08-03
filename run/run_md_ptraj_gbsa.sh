lig=$1
repeat=${2:-""} # optional parameter to run repeats, e.g. _1 _2 etc

echo "lig = " $lig "repeat = " $repeat

jobid=$(sbatch --parsable ~/scripts/1gpu.sh ~/cdk_scripts/standard_md.sh $lig "$repeat")

jobid=$(sbatch --dependency=afterok:$jobid --parsable ~/scripts/1cpu.sh ~/cdk_scripts/ptraj.py -i $lig -r "$repeat")

jobid=$(sbatch --dependency=afterok:$jobid --parsable ~/scripts/1cpu.sh ~/cdk_scripts/gbsa.py -i $lig -r "$repeat")





