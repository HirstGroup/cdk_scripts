#!/usr/bin/env python

import argparse
import os
import textwrap

from run import run
from util.delete_md2 import delete_mds


def center_strip(parm, traj, ignore_errors=False, print_cmd=False, repeat='', time='equi'):
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

    parm_root = os.path.splitext(parm)[0]
    traj_root = os.path.splitext(traj)[0]

    string = textwrap.dedent(f'''\n
    trajin {traj} 1 last 
    autoimage
    trajout {traj_root}_cent.nc        
    ''')

    with open(f'center{repeat}_{time}.ptraj', 'w') as f:
        f.write(string)

    run(f'cpptraj {parm} center{repeat}_{time}.ptraj', ignore_errors=ignore_errors, print_cmd=print_cmd)

    string = textwrap.dedent(f'''\
    trajin {traj_root}_cent.nc 
    strip :DMS:CL3:WAT:NA:CL:ETA:Na+:Cl-
    rms
    trajout {traj_root}_cent_strip.nc
    ''')

    with open(f'strip{repeat}_{time}.ptraj', 'w') as f:
        f.write(string)

    run(f'cpptraj {parm} strip{repeat}_{time}.ptraj', ignore_errors=ignore_errors, print_cmd=print_cmd)

    string = textwrap.dedent(f'''\
    parmstrip :DMS:CL3:WAT:NA:CL:ETA:Na+:Cl-
    parmwrite out {parm_root}_strip.parm7
    ''')

    with open(f'parmstrip{repeat}_{time}.ptraj', 'w') as f:
        f.write(string)

    run(f'cpptraj {parm} parmstrip{repeat}_{time}.ptraj', ignore_errors=ignore_errors, print_cmd=print_cmd)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run gbsa.py')

    # Required arguments
    parser.add_argument('--parm', help='Name of parameter file', required=True)
    parser.add_argument('--traj', help='Name of trajectory file', required=True)

    # Optional arguments 
    parser.add_argument('--cd', help='Name of directory to change into to run commands', required=False)
    parser.add_argument('--ignore_errors', help='Ignore errors when running system commands, i.e. continue and not crash, options YES/NO', required=False)
    parser.add_argument('-p', '--part', default='', help='Part pattern, e.g. 2, 3', required=False)
    parser.add_argument('-r', '--repeat', default='', help='Repeat pattern, e.g. _2, _3', required=False)
    parser.add_argument('-t', '--time', default='equi', help='Time pattern to run second MD, etc, e.g. equi, equi2', required=False)
    parser.add_argument('--test', default='NO', help='Test run', required=False)
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

    if args.test == 'YES':
        print_cmd = True
    else:
        print_cmd = False

    center_strip(args.parm, args.traj, ignore_errors=ignore_errors, print_cmd=print_cmd, repeat=args.repeat, time=f'{args.time}{args.part}')

    if args.delete:
        try:
            delete_mds(args.parm, args.traj, repeat=args.repeat, time=f'{args.time}{args.part}')
        except:
            if not ignore_errors:
                raise Exception('Errors occured in delete_mds, exit')
            else:
                print('Errors occured in delete_mds, continue')

    if args.cd:
        os.chdir('../')