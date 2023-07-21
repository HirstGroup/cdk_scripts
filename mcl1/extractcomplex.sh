#!/usr/bin/bash

# Extract protein pdb and ligand mol2 from complex prmtop and inpcrd files
# run from mcl1_inputs folder

set -e

lig1=$1
lig2=$2

cd L${lig1}/${lig1}-${lig2}_MD_NVT_rerun/

cat > extract-prot.ptraj <<EOF
trajin ${lig1}-${lig2}_merged.inpcrd
strip !:1-150
trajout ../l${lig1}-prot.pdb
EOF

cpptraj ${lig1}-${lig2}_merged.prmtop extract-prot.ptraj

cat > extract-lig.ptraj <<EOF
trajin ${lig1}-${lig2}_merged.inpcrd
strip !:L${lig1}
trajout ../l${lig1}.mol2
EOF

cpptraj ${lig1}-${lig2}_merged.prmtop extract-lig.ptraj

cd ../../




