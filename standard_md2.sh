#!/bin/bash

set -e

# default values
part=""
test=NO

while [[ $# > 0 ]]
do
key="$1"

case $key in
    -c|--complex)
    complex="$2"
    shift
    ;;
    -p|--part)
    part="$2"
    shift
    ;;
    -t|--time)
    time="$2"
    shift
    ;;
    --test)
    test=YES
    ;;
    *)
    echo "Unknown argument: $1"
    exit 1
    ;;
esac
shift
done

nstlim=5000000

if [ $test -eq "YES" ]
then
echo "Running test calculation"
nstlim=20000
fi


IFS='_' read -r receptor lig <<< "$complex"

cat > 08_equi.in << EOF
NPT MD w/No position restraints and PME (sander)
 &cntrl
  ntx    = 5,
  irest  = 1,
  ntpr   = 10000,
  ntwx   = 10000,
  ntwe   = 10000,
  ntwr   = 10000,
  ig     = -1,

  ntf    = 1,
  ntb    = 2,
  ntp = 1, pres0 = 1.0, taup = 10.0, gamma_ln = 1.0,
  cut    = 12.0,
  iwrap  = 1,
  nsnb   = 10,

  nstlim = $nstlim,
  t      = 0.0,
  nscm   = 1000,
  dt     = 0.002,

  temp0  = 300.0,
  tempi  = 300.0,
  ntt    = 3,
  tautp  = 2.0,

  ntc = 2, !Flag for SHAKE to perform bond length constraints. = 2 bonds involving hydrogen are constrained
  iwrap=1, ioutfm=1, ntwv=-1,ntave=1000,
&end
 &ewald
   skinnb=2, nfft1=96, nfft2=96, nfft3=96,
 / 
EOF

pmemd.cuda -O -i 08_equi.in -c ${complex}${repeat}_equi.rst7 -p ${complex}.parm7 -o ${complex}${repeat}_equi2.out -r ${complex}${repeat}_equi2.rst7 -x ${complex}${repeat}_equi2.nc -l ${complex}${repeat}_equi2.log