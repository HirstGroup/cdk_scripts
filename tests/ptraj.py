import argparse
import textwrap

from md import run


def center_strip(complex):
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

	string = textwrap.dedent('''\n
	trajin {complex}_{time}.nc 1 last 
	autoimage
	trajout {complex}_{time}_cent.nc		
	''')

	with open('center.ptraj', 'w') as f:
		f.write(string)

	run(f'cpptraj {complex}.parm7 center.ptraj')

	string = textwrap.dedent('''\
	trajin {complex}_{time}_cent.nc 
	strip :DMS:CL3:WAT:NA:CL:ETA:Na+:Cl-
	rms
	trajout {complex}_{time}_cent_strip.nc
	''')

	with open('strip.ptraj', 'w') as f:
		f.write(string)

	run(f'cpptraj {complex}.parm7 strip.ptraj')

	string = textwrap.dedent('''\
	parmstrip :DMS:CL3:WAT:NA:CL:ETA:Na+:Cl-
	parmwrite out {complex}_strip.parm7
	EOF
	''')

	with open('parmstrip.ptraj', 'w') as f:
		f.write(string)

	run(f'cpptraj {complex}.parm7 parmstrip.ptraj')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run gbsa.py')

    # Required arguments
    parser.add_argument('-i','--input', nargs='+', help='Name of complexes to run functions on',required=True)

    args = parser.parse_args()

    center_strip(complex)