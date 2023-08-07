set -e

sbatch=${1:-NO}

test1() {
# test for part 2

cp input/l23.parm7 output/
cp input/l23_equi.rst7 output/

cd output

bash -x ~/cdk_scripts/run/run_md_ptraj_gbsa.sh -c l23 -p 2 --sbatch $sbatch --test YES

ls gbsa2/l23_gbsa2.dat

cd ..
}

test2() {
# test for part 2 and repeat 2

cp input/l23.parm7 output/
cp input/l23_equi.rst7 output/l23_2_equi.rst7

cd output

bash -x ~/cdk_scripts/run/run_md_ptraj_gbsa.sh -c l23 -p 2 -r _2 --sbatch $sbatch --test YES

ls gbsa2_2/l23_gbsa2_2.dat

cd ..
}

test1
test2