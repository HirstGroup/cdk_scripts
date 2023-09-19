#!/bin/bash

set -e

# default values
part=NO
repeat=NO
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
    -r|--repeat)
    repeat="$2"
    shift
    ;;
    -t|--time)
    time="$2"
    shift
    ;;
    --test)
    test="$2"
    shift
    ;;
    *)
    echo "Unknown argument: $1"
    exit 1
    ;;
esac
shift
done

nstlim=5000000
ntpr=10000
if [ $test = "YES" ]; then
    echo "Running test calculation"
    ntpr=1000
    nstlim=2000
fi

if [ $repeat = "NO" ]; then
    repeat=""
fi 

if [ $part = "NO" ]; then
    part=""
    previous="${complex}${repeat}_heat2.rst7"
elif [ $part = 2 ]; then
    previous="${complex}${repeat}_equi.rst7"
else
    echo "Part greater different to NO or 2 not implemented"
    exit 1
fi

IFS='_' read -r receptor lig <<< "$complex"

cat > 08_equi.in << EOF
NPT MD w/No position restraints and PME (sander)
 &cntrl
  ntx    = 5,
  irest  = 1,
  ntpr   = $ntpr,
  ntwx   = $ntpr,
  ntwe   = $ntpr,
  ntwr   = $ntpr,
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

pmemd.cuda -O -i 08_equi.in -c $previous -p ${complex}.parm7 -o ${complex}${repeat}_equi${part}.out -r ${complex}${repeat}_equi${part}.rst7 -x ${complex}${repeat}_equi${part}.nc -l ${complex}${repeat}_equi${part}.log