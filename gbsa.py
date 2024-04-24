#!/usr/bin/env python

import argparse
import os
import textwrap

from run import run


def antegbsa(complex, ligandmask=None, part='', print_cmd=False, repeat=''):
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

    if ligandmask is None:
        if '_' in complex:
            ligandmask = complex.split('_')[1].upper()
        else:
            ligandmask = complex.upper()

    if not os.path.exists(f'gbsa{part}{repeat}'):
        os.makedirs(f'gbsa{part}{repeat}')

    os.chdir(f'gbsa{part}{repeat}')
    os.system(f'rm -f {complex}-complex.parm7 {complex}-receptor.parm7 {complex}-ligand.parm7')

    run(f'ante-MMPBSA.py -p ../{complex}.parm7 -c {complex}-complex.parm7 -r {complex}-receptor.parm7 -l {complex}-ligand.parm7 --strip-mask=:WAT:CL:NA --ligand-mask=:{ligandmask} --radii=mbondi2', print_cmd=print_cmd)

    os.chdir('../')


def gbsa(complex, part='', print_cmd=False, repeat=''):
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

    run(f'MMPBSA.py -O -i gbsa.in -o {complex}_gbsa{part}{repeat}.dat -cp {complex}-complex.parm7 -rp {complex}-receptor.parm7 -lp {complex}-ligand.parm7 -y ../{complex}{repeat}_equi{part}_cent_strip.nc', print_cmd=print_cmd)

    os.system('rm *.nc.0 reference.frc')    

    os.chdir('../')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run gbsa.py')

    # Required arguments
    parser.add_argument('-i','--input', help='Name of complexes to run functions on', required=True)

    # Optional arguments
    parser.add_argument('--cd', help='Name of directory to change into to run commands', required=False)
    parser.add_argument('-f','--functions', nargs='+', help='Function names to run', required=False)
    parser.add_argument('-l','--ligandmask', help='Ligandmask', required=False)
    parser.add_argument('-p', '--part', default='', help='Part pattern to run second MD, etc, e.g. 2, 3', required=False)
    parser.add_argument('-r','--repeat', default='', help='Repeat pattern, e.g. _2, _3', required=False)
    parser.add_argument('--test', default='NO', help='Test run', required=False)

    args = parser.parse_args()

    if args.cd:
        os.chdir(args.cd)

    complex = args.input

    if args.part == 'NO':
        args.part = ''

    if args.repeat == 'NO':
        args.repeat = ''

    if args.test == 'YES':
        print_cmd = True
    else:
        print_cmd = False

    if args.functions is None:
        antegbsa(complex, ligandmask=args.ligandmask, part=args.part, print_cmd=print_cmd, repeat=args.repeat)
        gbsa(complex, part=args.part, print_cmd=print_cmd, repeat=args.repeat)

    else:
        if 'antegbsa' in args.functions:
            antegbsa(complex, ligandmask=args.ligandmask, part=args.part, print_cmd=print_cmd, repeat=args.repeat)
        if 'gbsa' in args.functions:
            gbsa(complex, part=args.part, print_cmd=print_cmd, repeat=args.repeat)

    if args.cd:
        os.chdir('../')
