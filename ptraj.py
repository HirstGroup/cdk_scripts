#!/usr/bin/env python

import argparse
import os
import textwrap

from run import run


def center_strip(complex, time='equi', repeat=''):
	"""
	Center and strip MD trajectory

	Parameters
	----------
	complex: str
		name of complex

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
    parser.add_argument('-i','--input', help='Name of complex to run functions on',required=True)

    # Optional arguments
    parser.add_argument('--cd', help='Name of directory to change into to run commands',required=False)
    parser.add_argument('-r', '--repeat', default='', help='Repeat pattern, e.g. _2, _3', required=False)

    args = parser.parse_args()

    if args.cd:
    	os.chdir(args.cd)

    center_strip(args.input, repeat=args.repeat)

    if args.cd:
    	os.chdir('../')