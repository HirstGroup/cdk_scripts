#!/bib/bash

set -e

rm -rf output/*

cp input/l23.parm7 output/
cp input/l23_equi_cent.nc output/

cd output/

bash ../../run/run_ptraj_gbsa.sh -c l23 --ignore_errors YES

cd ..

