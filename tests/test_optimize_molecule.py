import filecmp
import os
import sys

sys.path.append('../')
from optimize_molecule import *

def compare_files(file1, file2, ignore_first):

	with open(file1) as f1, open(file2) as f2:

		lines1 = f1.readlines()[ignore_first:]
		lines2 = f2.readlines()[ignore_first:]

		return lines1 == lines2


def test_antechamber():

	os.chdir('output')

	os.system('pwd')

	run('cp ../input/1-1.mol2 .')

	antechamber('1-1.mol2')

	os.chdir('..')

	assert filecmp.cmp('input/1-1_gas.mol2', 'output/1-1_gas.mol2') is True

	assert filecmp.cmp('input/1-1_gas.frcmod', 'output/1-1_gas.frcmod') is True


def test_tleap():

	os.chdir('output')

	os.system('cp ../input/1-1_gas.mol2 .')
	os.system('cp ../input/1-1_gas.frcmod .')

	tleap('1-1_gas')

	os.chdir('..')

	assert compare_files('input/1-1_gas.parm7', 'output/1-1_gas.parm7', 1) is True
	assert compare_files('input/1-1_gas.rst7', 'output/1-1_gas.rst7', 1) is True	


def test_minimize():

	os.chdir('output')
	os.system('cp ../input/1-1_gas.parm7 .')
	os.system('cp ../input/1-1_gas.rst7 .')

	minimize('1-1_gas', '@1')

	os.chdir('..')


def test_make_bellymask_nonhydrogens():

	bellymask = make_bellymask_nonhydrogens('input/1-1_gas.mol2')

	assert bellymask == '@40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68'

test_minimize()

