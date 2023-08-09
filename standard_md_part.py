#!/usr/bin/env python

import argparse
import os
import textwrap

from run import run

def standard_md_part(complex, part, cd='NO', repeat='', test='NO'):

    if cd != 'NO':
        os.chdir(cd)

    if repeat == 'NO':
        repeat = ''

    nstlim = 5000000
    ntpr = 10000

    if test == 'YES':
        ntpr = 1000
        nstlim = 2000

    string = textwrap.dedent(f'''\
    NPT MD w/No position restraints and PME (sander)
     &cntrl
      ntx    = 5,
      irest  = 1,
      ntpr   = {ntpr},
      ntwx   = {ntpr},
      ntwe   = {ntpr},
      ntwr   = {ntpr},
      ig     = -1,

      ntf    = 1,
      ntb    = 2,
      ntp = 1, pres0 = 1.0, taup = 10.0, gamma_ln = 1.0,
      cut    = 12.0,
      iwrap  = 1,
      nsnb   = 10,

      nstlim = {nstlim},
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
    ''')

    with open(f'09_equi{part}{repeat}.in', 'w') as f:
        f.write(string)

    if part == 2:
        previous_rst7 = f'{complex}{repeat}_heat2.rst7'
    else:
        previous_rst7 = f'{complex}{repeat}_equi{part-1}.rst7'

    run(f'pmemd.cuda -O -i 09_equi{part}{repeat}.in -c {previous_rst7} -p {complex}.parm7 -o {complex}{repeat}_equi{part}.out -r {complex}{repeat}_equi{part}.rst7 -x {complex}{repeat}_equi{part}.nc -l {complex}{repeat}_equi{part}.log')

    if cd:
        os.chdir('../')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run a part of an MD (i.e. sequential MDs)')

    # Required arguments
    parser.add_argument('-c','--complex', help='Name of complexes to run', required=True)

    # Optional arguments
    parser.add_argument('--cd', default='NO', help='Name of directory to change into to run commands', required=False)
    parser.add_argument('-p', '--part', type=int, help='Part number to run MD, etc, e.g. 2, 3', required=False)
    parser.add_argument('-r','--repeat', default='NO', help='Repeat pattern, e.g. _2, _3', required=False)
    parser.add_argument('--test', default='NO', help='Do test run, i.e. very short MD', required=False)

    args = parser.parse_args()

    standard_md_part(complex=args.complex, part=args.part, cd=args.cd, repeat=args.repeat, test=args.test)