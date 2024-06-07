import filecmp
import sys
import textwrap

sys.path.append('../')
from gbsa_decomp import *

from tests.compare_files import compare_files


def test_antegbsa():

	os.system('rm -rf output/*')

	os.system('cp input/a01.parm7 output/')

	os.chdir('output')

	antegbsa('a01')

	os.chdir('../')

	for file in ['a01-complex.parm7', 'a01-ligand.parm7', 'a01-receptor.parm7']:

		compare_files(f'input/{file}', f'output/gbsa_decomp/{file}', 1)


def test_gbsa_decomp():

	os.system('rm -rf output/*')

	os.system('cp input/a01.parm7 output/')
	os.system('cp input/a01_equi_cent_strip.nc output/')

	os.chdir('output')

	antegbsa('a01')
	gbsa_decomp('a01')

	os.chdir('../')

	compare_files('input/a01_gbsa_decomp.dat', 'output/gbsa_decomp/a01_gbsa_decomp.dat', 1)