set -e

cp input/l23.parm7 output/
cp input/l23_equi.rst7 output/

cd output

bash ~/cdk_scripts/run/run_md_ptraj_gbsa.sh -c l23 -p 2 --test

cd ..
