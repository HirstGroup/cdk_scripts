set -e

test1() {
# test for part 2

cp input/l23.parm7 output/
cp input/l23_equi.rst7 output/

cd output

bash -x ~/cdk_scripts/run/run_md_ptraj_gbsa_part2.sh -c l23 --test YES

cd ..
}

test2() {
# test for part 2 and repeat 2

cp input/l23.parm7 output/
cp input/l23_equi.rst7 output/l23_2_equi.rst7

cd output

bash -x ~/cdk_scripts/run/run_md_ptraj_gbsa_part2.sh -c l23 -r _2 --test YES

cd ..
}

test1
test2

ls output/gbsa2/l23_gbsa2.dat
ls output/gbsa2_2/l23_gbsa2_2.dat