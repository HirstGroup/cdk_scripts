import filecmp
import sys

sys.path.append('../')
from gbsa import *

def test_get_charge():

	os.system('ls input/1-1.mol2')

	charge = get_charge('input/1-1.mol2')

	assert charge == 2


def test_create_resp1_file():

	charge = get_charge('input/1-1.mol2')

	create_resp1_file('input/1-1.mol2', 'output/1-1_opt.gau', charge, cpu=10)

	assert filecmp.cmp('input/1-1_opt.gau', 'output/1-1_opt.gau') is True


def test_check_resp1_output():

	assert check_resp1_output('input/1-1_opt.log', 'PHXJVRSECIGDHY-UHFFFAOYSA-P') == 'OK'


def test_get_inchikey():

	assert get_inchikey('input/1-1.mol2') == 'PHXJVRSECIGDHY-UHFFFAOYSA-P'


def test_create_resp2_file():

	create_resp2_file('input/1-1_opt.log', 'output/1-1_esp.gau', charge=2, cpu=10)

	assert filecmp.cmp('input/1-1_esp.gau', 'output/1-1_esp.gau') is True


def test_create_resp3_file():

	create_resp3_file('input/1-1_esp.log', 'output/1-1_resp3.mol2', 'output/1-1_resp.mol2', 'input/1-1.mol2', 'LIG')

	assert filecmp.cmp('input/1-1_resp.mol2', 'output/1-1_resp.mol2') is True


test_check_resp3_file()