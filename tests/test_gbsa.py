import filecmp
import sys
import textwrap

sys.path.append('../')
from gbsa import *

from tests.compare_files import compare_files


def test_antegbsa():

	os.system('cp input/a01.parm7 output/')

	os.chdir('output')

	antegbsa('a01')

	os.chdir('../')

	for file in ['a01-complex.parm7', 'a01-ligand.parm7', 'a01-receptor.parm7']:

		compare_files(f'input/{file}', f'output/gbsa/{file}', 1)


def test_gbsa():

	os.system('cp input/a01.parm7 output/')
	os.system('cp input/a01_equi_cent_strip.nc output/')

	os.chdir('output')

	antegbsa('a01')
	gbsa('a01')

	os.chdir('../')

	compare_files('input/a01_gbsa.dat', 'output/gbsa/a01_gbsa.dat', 1)


def test_gbsa_arg():

	os.system('cp input/a01.parm7 output/')
	os.system('cp input/a01_equi_cent_strip.nc output/')
	
	os.system('rm -rf output/gbsa/')

	os.chdir('output')

	os.system('python ../../gbsa.py -i a01')

	os.chdir('../')

	compare_files('input/a01_gbsa.dat', 'output/gbsa/a01_gbsa.dat', 1)


def test_gbsa_arg2():

	os.system('cp input/a01.parm7 output/')
	os.system('cp input/a01_equi_cent_strip.nc output/')

	os.system('rm -rf output/gbsa/')

	os.system('python ../gbsa.py -i a01 --cd output')

	compare_files('input/a01_gbsa.dat', 'output/gbsa/a01_gbsa.dat', 1)


def test_gbsa_main():
	# test gbsa main with part option

    os.system('cp input/l23.parm7 output/')
    os.system('cp input/l23_equi2_cent_strip.nc output/')

    os.system('rm -rf output/gbsa2/')

    os.system('python ../gbsa.py -i l23 --cd output -p 2')

    compare_files('input/l23_gbsa2.dat', 'output/gbsa2/l23_gbsa2.dat', 1)

