#!/bin/bash

set -e

complex=$1
repeat=$2 # optional parameter to run repeats, e.g. _1 _2 etc

IFS='_' read -r receptor lig <<< "$complex"

cat > 01_min.in << EOF
#NVT MD w/No position restraints and PME (sander)

 &cntrl
  imin   = 1,
  maxcyc = 20000,
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
  
  ntc=2,
  
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
  maxcyc = 20000,
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
  
  ntc = 2,
  
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
  maxcyc = 20000,
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
  
  ntc = 2,
  
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
  maxcyc = 20000,
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
  
  ntc = 2,
  
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
  imin   = 1,
  maxcyc = 20000,
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
  
  ntc = 2,
  
 &end
  &ewald
  skinnb = 2,
  nfft1 = 96, nfft2 = 96, nfft3 = 96,
 /  
EOF

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

  nstlim = 500000, 
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

  nstlim = 500000, 
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

cat > 08_equi.in << EOF
NPT MD w/No position restraints and PME (sander)
 &cntrl
  ntx    = 5,
  irest  = 1,
  ntpr   = 1000,
  ntwx   = 1000,
  ntwe   = 1000,
  ntwr   = 1000,
  ig     = -1,

  ntf    = 1,
  ntb    = 2,
  ntp = 1, pres0 = 1.0, taup = 10.0, gamma_ln = 1.0,
  cut    = 12.0,
  iwrap  = 1,
  nsnb   = 10,

  nstlim = 5000000,
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