#!/usr/bin/env python

import argparse
import os
import textwrap

from util.delete_md import delete_mds
from run import run


def center_strip(complex, ignore_errors=False, print_cmd=False, repeat='', time='equi'):
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

    with open(f'center{repeat}_{time}.ptraj', 'w') as f:
        f.write(string)

    run(f'cpptraj {complex}.parm7 center{repeat}_{time}.ptraj', ignore_errors=ignore_errors, print_cmd=print_cmd)

    string = textwrap.dedent(f'''\
    trajin {complex}{repeat}_{time}_cent.nc 
    strip :DMS:CL3:WAT:NA:CL:ETA:Na+:Cl-
    rms
    trajout {complex}{repeat}_{time}_cent_strip.nc
    ''')

    with open(f'strip{repeat}_{time}.ptraj', 'w') as f:
        f.write(string)

    run(f'cpptraj {complex}.parm7 strip{repeat}_{time}.ptraj', ignore_errors=ignore_errors, print_cmd=print_cmd)

    string = textwrap.dedent(f'''\
    parmstrip :DMS:CL3:WAT:NA:CL:ETA:Na+:Cl-
    parmwrite out {complex}_strip.parm7
    ''')

    with open(f'parmstrip{repeat}_{time}.ptraj', 'w') as f:
        f.write(string)

    run(f'cpptraj {complex}.parm7 parmstrip{repeat}_{time}.ptraj', ignore_errors=ignore_errors, print_cmd=print_cmd)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run gbsa.py')

    # Required arguments
    parser.add_argument('-i','--input', help='Name of complex to run functions on', required=True)

    # Optional arguments 
    parser.add_argument('--cd', help='Name of directory to change into to run commands', required=False)
    parser.add_argument('--ignore_errors', help='Ignore errors when running system commands, i.e. continue and not crash, options YES/NO', required=False)
    parser.add_argument('-p', '--part', default='', help='Part pattern, e.g. 2, 3', required=False)
    parser.add_argument('-r', '--repeat', default='', help='Repeat pattern, e.g. _2, _3', required=False)
    parser.add_argument('-t', '--time', default='equi', help='Time pattern to run second MD, etc, e.g. equi, equi2', required=False)
    parser.add_argument('--delete', action='store_true', help='Delete original MD trajectory files', required=False)

    args = parser.parse_args()

    if args.cd:
        os.chdir(args.cd)

    if args.ignore_errors == 'YES':
        ignore_errors = True
    else:
        ignore_errors = False

    if args.part == 'NO':
        args.part = ''

    if args.repeat == 'NO':
        args.repeat = ''

    center_strip(args.input, ignore_errors=ignore_errors, repeat=args.repeat, time=f'{args.time}{args.part}')

    if args.delete:
        try:
            delete_mds(args.input, repeat=args.repeat, time=f'{args.time}{args.part}')
        except:
            if not ignore_errors:
                raise Exception('Errors occured in delete_mds, exit')
            else:
                print('Errors occured in delete_mds, continue')

    if args.cd:
        os.chdir('../')
