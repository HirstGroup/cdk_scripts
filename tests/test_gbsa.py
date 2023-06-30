import filecmp
import sys
import textwrap

sys.path.append('../')
from gbsa import *

from compare_files import compare_files

def test_antegbsa():

	os.system('cp input/a01.parm7 output/')

	os.chdir('output')

	antegbsa('a01')

	os.chdir('../')

	for file in ['a01-complex.parm7', 'a01-ligand.parm7', 'a01-receptor.parm7']:

		compare_files(f'input/{file}', f'output/gbsa/{file}', 1)


def test_gbsa():

	os.chdir('output')

	gbsa('a01')

	os.chdir('../')

	compare_files('input/a01_gbsa.dat', 'output/gbsa/a01_gbsa.dat', 1)