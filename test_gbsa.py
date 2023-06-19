import filecmp

from gbsa import *

def test_get_charge():

	os.system('ls tests/1-1.mol2')

	charge = get_charge('tests/1-1.mol2')

	assert charge == 2

def test_create_resp1_file():

	charge = get_charge('tests/1-1.mol2')

	create_resp1_file('tests/1-1.mol2', 'tests/output/1-1_opt.gau', charge, cpu=10)

	assert filecmp.cmp('tests/1-1_opt.gau', 'tests/output/1-1_opt.gau') is True


test_create_resp1_file()