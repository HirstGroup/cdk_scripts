#!/bin/bash

# shorten an MD in intervals of 10 and rename original

complex=$1
repeat=${2:-""} # optional parameter to run repeats, e.g. _1 _2 etc

mv ${complex}${repeat}_equi_cent_strip.nc ${complex}${repeat}_equi_cent_strip_original.nc

cat > shorten.ptraj << EOF
trajin ${complex}${repeat}_equi_cent_strip_original.nc 1 last 10
trajout ${complex}${repeat}_equi_cent_strip.nc
EOF

cpptraj ${complex}_strip.parm7 shorten.ptraj
