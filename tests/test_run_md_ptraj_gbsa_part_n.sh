set -e

rm -rf output/*
cp input/l23.parm7 output/
cp input/l23_heat2.rst7 output/l23_3_equi.rst7

cd input/

python ../../run/run_md_ptraj_gbsa_part_n.py -c l23 -p 2 3 4 5 6 7 8 9 10 -r _3 --test YES

cd ../