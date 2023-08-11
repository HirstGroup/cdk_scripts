set -e

# TEST 1

first_part=2
last_part=10

rm -rf output/*
cp input/l23.parm7 output/
cp input/l23_heat2.rst7 output/l23_3_equi.rst7

cd output/

python ../../run/run_md_ptraj_gbsa_part_n.py -c l23 -fp $first_part -lp $last_part -r _3 --test YES

for i in $(seq $first_part $last_part); do
ls l23_3_equi$i.out l23_3_equi$i.rst7 l23_3_equi${i}_cent_strip.nc gbsa${i}_3/l23_gbsa${i}_3.dat
done

<<<<<<< HEAD
cd ../
=======
python ../../run/run_md_ptraj_gbsa_part_n.py -c l23 -p 2 3 4 5 6 7 8 9 10 -r _3 --test YES
>>>>>>> 3caf28c7927efb3efc9ed31cde0e362dbae74419

echo "TEST 1 PASSED"