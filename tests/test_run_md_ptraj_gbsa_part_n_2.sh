set -e

# TEST 2 (with sequential option)

first_part=2
last_part=3

rm -rf output/*
cp input/l23.parm7 output/
cp input/l23_heat2.rst7 output/l23_3_equi.rst7

cd output/

mkdir test

python ../../run/run_md_ptraj_gbsa_part_n.py -c l23 -fp $first_part -lp $last_part -r _3 --sequential --test YES

for i in $(seq $first_part $last_part); do
ls l23_3_equi$i.out l23_3_equi$i.rst7 l23_3_equi${i}_cent_strip.nc gbsa${i}_3/l23_gbsa${i}_3.dat
mv l23_3_equi$i.out l23_3_equi$i.rst7 l23_3_equi${i}_cent_strip.nc gbsa${i}_3/l23_gbsa${i}_3.dat test/
done

cd ..

echo "TEST 2 PASSED"