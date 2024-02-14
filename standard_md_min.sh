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
  
  ntc = 2,
  
 &end
  &ewald
  skinnb = 2,
  nfft1 = 96, nfft2 = 96, nfft3 = 96,
 /  
EOF


pmemd.cuda -O -i 01_min.in -c ${complex}.rst7 -p ${complex}.parm7 -o ${complex}${repeat}_min1.out -r ${complex}${repeat}_min1.rst7 -x ${complex}${repeat}_min1.nc -ref ${complex}.rst7
pmemd.cuda -O -i 02_min.in -c ${complex}${repeat}_min1.rst7 -p ${complex}.parm7 -o ${complex}${repeat}_min2.out -r ${complex}${repeat}_min2.rst7 -x ${complex}${repeat}_min2.nc -ref ${complex}${repeat}_min1.rst7
pmemd.cuda -O -i 03_min.in -c ${complex}${repeat}_min2.rst7 -p ${complex}.parm7 -o ${complex}${repeat}_min3.out -r ${complex}${repeat}_min3.rst7 -x ${complex}${repeat}_min3.nc -ref ${complex}${repeat}_min2.rst7
pmemd.cuda -O -i 04_min.in -c ${complex}${repeat}_min3.rst7 -p ${complex}.parm7 -o ${complex}${repeat}_min4.out -r ${complex}${repeat}_min4.rst7 -x ${complex}${repeat}_min4.nc -ref ${complex}${repeat}_min3.rst7
pmemd.cuda -O -i 05_min.in -c ${complex}${repeat}_min4.rst7 -p ${complex}.parm7 -o ${complex}${repeat}_min5.out -r ${complex}${repeat}_min5.rst7 -x ${complex}${repeat}_min5.nc -ref ${complex}${repeat}_min4.rst7