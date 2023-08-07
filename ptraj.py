#!/usr/bin/env python

import argparse
import os
import textwrap

from util.delete_md import delete_mds
from run import run


def center_strip(complex, repeat='', time='equi'):
    """
    Center and strip MD trajectory

    Parameters
    ----------
    complex: str
        name of complex
    part : srt
        part pattern, e.g. 2, 3
    time : str
        time pattern, e.g. equi, equi2

    Returns
    -------
    None (creates output trajectory files)
    """

    string = textwrap.dedent(f'''\n
    trajin {complex}{repeat}_{time}.nc 1 last 
    autoimage
    trajout {complex}{repeat}_{time}_cent.nc        
    ''')

    with open('center.ptraj', 'w') as f:
        f.write(string)

    run(f'cpptraj {complex}.parm7 center.ptraj')

    string = textwrap.dedent(f'''\
    trajin {complex}{repeat}_{time}_cent.nc 
    strip :DMS:CL3:WAT:NA:CL:ETA:Na+:Cl-
    rms
    trajout {complex}{repeat}_{time}_cent_strip.nc
    ''')

    with open('strip.ptraj', 'w') as f:
        f.write(string)

    run(f'cpptraj {complex}.parm7 strip.ptraj')

    string = textwrap.dedent(f'''\
    parmstrip :DMS:CL3:WAT:NA:CL:ETA:Na+:Cl-
    parmwrite out {complex}_strip.parm7
    ''')

    with open('parmstrip.ptraj', 'w') as f:
        f.write(string)

    run(f'cpptraj {complex}.parm7 parmstrip.ptraj')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run gbsa.py')

    # Required arguments
    parser.add_argument('-i','--input', help='Name of complex to run functions on', required=True)

    # Optional arguments 
    parser.add_argument('--cd', help='Name of directory to change into to run commands', required=False)
    parser.add_argument('-p', '--part', default='', help='Part pattern, e.g. 2, 3', required=False)
    parser.add_argument('-r', '--repeat', default='', help='Repeat pattern, e.g. _2, _3', required=False)
    parser.add_argument('-t', '--time', default='equi', help='Time pattern to run second MD, etc, e.g. equi, equi2', required=False)
    parser.add_argument('--delete', action='store_true', help='Delete original MD trajectory files', required=False)

    args = parser.parse_args()

    if args.cd:
        os.chdir(args.cd)

    if args.part == 'NO':
        args.part = ''

    if args.repeat == 'NO':
        args.repeat = ''

    center_strip(args.input, repeat=args.repeat, time=args.time + args.part)

    if args.delete:
        delete_mds(args.input, repeat=args.repeat, time=args.time + args.part)

    if args.cd:
        os.chdir('../')