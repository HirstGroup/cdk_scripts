import argparse
import os
import textwrap

from md import run

def antegbsa(complex):
	"""
	Run ante-MMPBSA.py to prepare topolopy files for MMPBSA.py

	Parameters
	----------
	complex: str
		complex name

	Returns
	-------
	None (runs ante-MMPBSA.py)
	"""

	ligandmask = complex.upper()

	if not os.path.exists('gbsa'):
    	os.makedirs('gbsa')

	os.chdir('gbsa')
	os.system(f'rm -f {complex}-complex.parm7 {complex}-receptor.parm7 {complex}-ligand.parm7')

	run(f'ante-MMPBSA.py -p {complex}.parm7 -c {complex}-complex.parm7 -r {complex}-receptor.parm7 -l {complex}-ligand.parm7 --strip-mask=:WAT:CL:NA --ligand-mask={ligandmask} --radii=mbondi2')

	os.chdir('../')

def gbsa(complex):
	"""
	Run MMPBSA.py

	Parameters
	----------
	complex: str
		complex name

	Returns
	-------
	None (runs MMPBSA.py)
	"""

	os.chdir('gbsa')

	string = textwrap.dedent('''\
	&general 
	  interval=1, keep_files=1, netcdf=1,
	/   
	&gb 
	  igb=5, saltcon=0.100,
	/
	''')

	with open('gbsa.in', 'w') as f:
		f.write(string)

	run(f'MMPBSA.py -O -i gbsa.in -o {complex}_gbsa.dat -cp {complex}-complex.parm7 -rp {complex}-receptor.parm7 -lp {complex}-ligand.parm7 -y {complex}_equi_cent_strip.nc')

	os.chdir('../')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run gbsa.py')

    # Required arguments
    parser.add_argument('-f','--function', help='Function name to run',required=True)
    parser.add_argument('-i','--input', nargs='+', help='Name of complexes to run functions on',required=True)

    args = parser.parse_args()

    if args.function == 'antegbsa':
    	antegbsa(complex)
    if args.function == 'gbsa':
    	gbsa(complex)
