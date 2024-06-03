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

if [ $repeat = "NO" ]; then
repeat=""
fi 

maxcyc=20000

if [ $test = "YES" ]; then
maxcyc=100
fi

cat > 01_min.in << EOF
#NVT MD w/No position restraints and PME (sander)

 &cntrl
  imin   = 1,
  maxcyc = $maxcyc,
  ntmin  = 2,
  ntpr   = 1000,
  ntwx   = 1000,
  ntwe   = 1000,
  ntwr   = 1000,

  cut    = 12.0,
  iwrap  = 1,
  nsnb   = 10,

  iwrap  = 1,
  ioutfm = 1,
  
  ntr=1,
  restraintmask = '!(:WAT,NA,CL)',
  restraint_wt = 100.0

 &end
  &ewald
  skinnb = 2,
  nfft1 = 96, nfft2 = 96, nfft3 = 96,
 /  
EOF

cat > 02_min.in << EOF
#NVT MD w/No position restraints and PME (sander)

 &cntrl
  imin   = 1,
  maxcyc = $maxcyc,
  ntmin  = 2,
  ntpr   = 1000,
  ntwx   = 1000,
  ntwe   = 1000,
  ntwr   = 1000,

  cut    = 12.0,
  iwrap  = 1,
  nsnb   = 10,

  iwrap  = 1,
  ioutfm = 1,
  
  ntr = 1,
  restraintmask = '!(:WAT,NA,CL) & !@H=',
  restraint_wt = 100.0

 &end
  &ewald
  skinnb = 2,
  nfft1 = 96, nfft2 = 96, nfft3 = 96,
 /  
EOF

cat > 03_min.in << EOF
#NVT MD w/No position restraints and PME (sander)

 &cntrl
  imin   = 1,
  maxcyc = $maxcyc,
  ntmin  = 2,
  ntpr   = 1000,
  ntwx   = 1000,
  ntwe   = 1000,
  ntwr   = 1000,

  cut    = 12.0,
  iwrap  = 1,
  nsnb   = 10,

  iwrap  = 1,
  ioutfm = 1,
  
  ntr = 1,
  restraintmask = '@CA,C,O,N',
  restraint_wt = 100.0

 &end
  &ewald
  skinnb = 2,
  nfft1 = 96, nfft2 = 96, nfft3 = 96,
 /  
EOF

cat > 04_min.in << EOF
#NVT MD w/No position restraints and PME (sander)

 &cntrl
  imin   = 1,
  maxcyc = $maxcyc,
  ntmin  = 2,
  ntpr   = 1000,
  ntwx   = 1000,
  ntwe   = 1000,
  ntwr   = 1000,

  cut    = 12.0,
  iwrap  = 1,
  nsnb   = 10,

  iwrap  = 1,
  ioutfm = 1,
  
  ntr = 1,
  restraintmask = '@CA,C,O',
  restraint_wt = 100.0

 &end
  &ewald
  skinnb = 2,
  nfft1 = 96, nfft2 = 96, nfft3 = 96,
 /  
EOF

cat > 05_min.in << EOF
#NVT MD w/No position restraints and PME (sander)

 &cntrl
  imin   = 1, !Flag to run minimization, =1 ON
  maxcyc = $maxcyc,
  ntmin  = 1, !Flag for method of minimization, =1 Steepest descent for ncyc cycles, then conjugate gradient
  ncyc   = 1000, !minimization will be switched from steepest descent to conjugate gradient after NCYC cycles
  ntpr   = 1000,
  ntwx   = 1000,
  ntwe   = 1000,
  ntwr   = 1000,

  cut    = 12.0,
  iwrap  = 1,
  nsnb   = 10, !Determines the frequency of nonbonded list updates when igb=0 and nbflag=0, default is 25

  iwrap  = 1,
  ioutfm = 1,
  
 &end
  &ewald
  skinnb = 2, !Width of nonbonded "skin", default is 2 
  nfft1 = 96, nfft2 = 96, nfft3 = 96, !Size of the charge grid
 /  
EOF

nstlim=500000

if [ $test = "YES" ]; then
nstlim=100
fi

cat > 06_heat.in << EOF
NVT MD w/No position restraints and PME (sander)
 &cntrl
  ntx    = 1,
  irest  = 0,
  ntpr   = 10000,
  ntwx   = 10000,
  ntwe   = 10000,
  ntwr   = 10000,
  ig     = -1,

  ntf    = 2, !Force evaluation. = 2 bond interactions involving H-atoms omitted (use with NTC=2)
  ntb    = 1, != 1 constant volume
  cut    = 12.0,
  iwrap  = 1,
  nsnb   = 10, !determines the frequency of non-bonded list updates when igb=0 and nbflag=0

  nstlim = $nstlim, 
  t      = 0.0,
  nscm   = 1000, !Flag for the removal of translational and rotational center-of-mass (COM) motion at regular intervals (default is 1000)
  dt     = 0.002,

  temp0  = 300.0,
  tempi  = 0.0,
  ntt    = 3,
  tautp  = 2.0,
  gamma_ln = 2.0,

  ntr=1,
  restraintmask = '!(:WAT,NA,CL)',
  restraint_wt = 5.0,

  ntc = 2, !Flag for SHAKE to perform bond length constraints. = 2 bonds involving hydrogen are constrained
  iwrap = 1, ioutfm = 1, ntwv = -1, ntave = 1000,
&end
 &ewald
   skinnb = 2, nfft1 = 96, nfft2 = 96, nfft3 = 96,
 /
EOF

cat > 07_heat.in << EOF
NVT MD w/No position restraints and PME (sander)
 &cntrl
  ntx    = 1,
  irest  = 0,
  ntpr   = 10000,
  ntwx   = 10000,
  ntwe   = 10000,
  ntwr   = 10000,
  ig     = -1,

  ntf    = 2, !Force evaluation. = 2 bond interactions involving H-atoms omitted (use with NTC=2)
  ntb    = 1, != 1 constant volume
  cut    = 12.0,
  iwrap  = 1,
  nsnb   = 10, !determines the frequency of non-bonded list updates when igb=0 and nbflag=0

  nstlim = $nstlim, 
  t      = 0.0,
  nscm   = 1000, !Flag for the removal of translational and rotational center-of-mass (COM) motion at regular intervals (default is 1000)
  dt     = 0.002,

  temp0  = 300.0,
  tempi  = 0.0,
  ntt    = 3,
  tautp  = 2.0,
  gamma_ln = 2.0,

  ntc = 2, !Flag for SHAKE to perform bond length constraints. = 2 bonds involving hydrogen are constrained
  iwrap = 1, ioutfm = 1, ntwv = -1, ntave = 1000,
&end
 &ewald
   skinnb = 2, nfft1 = 96, nfft2 = 96, nfft3 = 96,
 /
EOF

ntwx=10000
nstlim=5000000

if [ $test = "YES" ]; then
ntwx=100
nstlim=100
fi

cat > 08_equi.in << EOF
NPT MD w/No position restraints and PME (sander)
 &cntrl
  ntx    = 5,
  irest  = 1,
  ntpr   = 10000,
  ntwx   = $ntwx,
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

pmemd.cuda -O -i 01_min.in -c ${complex}.rst7 -p ${complex}.parm7 -o ${complex}${repeat}_min1.out -r ${complex}${repeat}_min1.rst7 -x ${complex}${repeat}_min1.nc -ref ${complex}.rst7
pmemd.cuda -O -i 02_min.in -c ${complex}${repeat}_min1.rst7 -p ${complex}.parm7 -o ${complex}${repeat}_min2.out -r ${complex}${repeat}_min2.rst7 -x ${complex}${repeat}_min2.nc -ref ${complex}${repeat}_min1.rst7
pmemd.cuda -O -i 03_min.in -c ${complex}${repeat}_min2.rst7 -p ${complex}.parm7 -o ${complex}${repeat}_min3.out -r ${complex}${repeat}_min3.rst7 -x ${complex}${repeat}_min3.nc -ref ${complex}${repeat}_min2.rst7
pmemd.cuda -O -i 04_min.in -c ${complex}${repeat}_min3.rst7 -p ${complex}.parm7 -o ${complex}${repeat}_min4.out -r ${complex}${repeat}_min4.rst7 -x ${complex}${repeat}_min4.nc -ref ${complex}${repeat}_min3.rst7
pmemd.cuda -O -i 05_min.in -c ${complex}${repeat}_min4.rst7 -p ${complex}.parm7 -o ${complex}${repeat}_min5.out -r ${complex}${repeat}_min5.rst7 -x ${complex}${repeat}_min5.nc -ref ${complex}${repeat}_min4.rst7

pmemd.cuda -O -i 06_heat.in -c ${complex}${repeat}_min5.rst7 -p ${complex}.parm7 -o ${complex}${repeat}_heat1.out -r ${complex}${repeat}_heat1.rst7 -x ${complex}${repeat}_heat1.nc -ref ${complex}${repeat}_min5.rst7

pmemd.cuda -O -i 07_heat.in -c ${complex}${repeat}_heat1.rst7 -p ${complex}.parm7 -o ${complex}${repeat}_heat2.out -r ${complex}${repeat}_heat2.rst7 -x ${complex}${repeat}_heat2.nc -ref ${complex}${repeat}_heat1.rst7

pmemd.cuda -O -i 08_equi.in -c ${complex}${repeat}_heat2.rst7 -p ${complex}.parm7 -o ${complex}${repeat}_equi.out -r ${complex}${repeat}_equi.rst7 -x ${complex}${repeat}_equi.nc -l ${complex}${repeat}_equi.log