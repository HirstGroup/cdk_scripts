#!/usr/bin/env python

import argparse
import os
import textwrap

from run import run


def antegbsa(complex, part='', repeat=''):
    """
    Run ante-MMPBSA.py to prepare topolopy files for MMPBSA.py

    Parameters
    ----------
    complex: str
        complex name
    part : srt
        part pattern, e.g. 2, 3
    repeat : str
        repeat pattern, e.g. _2, _3

    Returns
    -------
    None (runs ante-MMPBSA.py)
    """

    ligandmask = complex.upper()

    if not os.path.exists(f'gbsa{part}{repeat}'):
        os.makedirs(f'gbsa{part}{repeat}')

    os.chdir(f'gbsa{part}{repeat}')
    os.system(f'rm -f {complex}-complex.parm7 {complex}-receptor.parm7 {complex}-ligand.parm7')

    run(f'ante-MMPBSA.py -p ../{complex}.parm7 -c {complex}-complex.parm7 -r {complex}-receptor.parm7 -l {complex}-ligand.parm7 --strip-mask=:WAT:CL:NA --ligand-mask=:{ligandmask} --radii=mbondi2')

    os.chdir('../')


def gbsa(complex, part='', repeat=''):
    """
    Run MMPBSA.py

    Parameters
    ----------
    complex: str
        complex name
    part : srt
        part pattern, e.g. 2, 3
    repeat : str
        repeat pattern, e.g. _2, _3

    Returns
    -------
    None (runs MMPBSA.py)
    """

    os.chdir(f'gbsa{part}{repeat}')

    string = textwrap.dedent('''\
    &general 
      interval=1, netcdf=1,
    /   
    &gb 
      igb=5, saltcon=0.100,
    /
    ''')

    with open('gbsa.in', 'w') as f:
        f.write(string)

    run(f'MMPBSA.py -O -i gbsa.in -o {complex}_gbsa{part}{repeat}.dat -cp {complex}-complex.parm7 -rp {complex}-receptor.parm7 -lp {complex}-ligand.parm7 -y ../{complex}{repeat}_equi{part}_cent_strip.nc')

    os.chdir('../')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run gbsa.py')

    # Required arguments
    parser.add_argument('-i','--input', help='Name of complexes to run functions on', required=True)

    # Optional arguments
    parser.add_argument('--cd', help='Name of directory to change into to run commands', required=False)
    parser.add_argument('-f','--functions', nargs='+', help='Function names to run', required=False)
    parser.add_argument('-p', '--part', default='', help='Part pattern to run second MD, etc, e.g. 2, 3', required=False)
    parser.add_argument('-r','--repeat', default='', help='Repeat pattern, e.g. _2, _3', required=False)

    args = parser.parse_args()

    if args.cd:
        os.chdir(args.cd)

    complex = args.input

    if args.functions is None:
        antegbsa(complex, part=args.part, repeat=args.repeat)
        gbsa(complex, part=args.part, repeat=args.repeat)

    else:
        if 'antegbsa' in args.functions:
            antegbsa(complex, args.repeat)
        if 'gbsa' in args.functions:
            gbsa(complex, args.part, args.repeat)

    if args.cd:
        os.chdir('../')
