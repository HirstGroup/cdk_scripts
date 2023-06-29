import filecmp
import sys
import textwrap

sys.path.append('../')
from md import *


def test_run():

	output = run('echo Hello')

	assert output == 'Hello\n'

	try:
		output = run('command_not_found')
	except:
		print('Exception handled')

	output = run('echo End')


def test_get_charge():

	os.system('ls input/1-1.mol2')

	charge = get_charge('input/1-1.mol2')

	assert charge == 2


def test_get_charge2():

	charge = get_charge('input/a01_confs_dock_best_0_p.mol2')

	print(charge)

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


def test_make_tleap_input1():

	ligandname = '../input/a01_resp'

	ligand = '../input/a01_resp.mol2'

	receptor = '../input/6td3_E_F.pdb'

	complex = '6td3_a01'

	neutral_expression = ''

	tleap_input = make_tleap_input(ligandname, ligand, receptor, complex, neutral_expression)

	output = textwrap.dedent(f'''\
	source leaprc.protein.ff19SB
	source leaprc.phosaa19SB
	source leaprc.gaff
	loadamberparams ../input/a01_resp.frcmod
	loadamberparams frcmod.ionsjc_tip3p
	source leaprc.water.tip3p
	LIG = loadMol2 ../input/a01_resp.mol2
	receptor = loadPDB ../input/6td3_E_F.pdb
	complex = combine {{receptor LIG}}
	set default PBRadii mbondi2

	solvatebox complex TIP3PBOX 10.0
	savepdb complex 6td3_a01-box.pdb
	saveAmberParm complex 6td3_a01.parm7 6td3_a01.rst7
	quit
	''')	
	assert tleap_input == output


def test_make_tleap_input2():

	os.chdir('output')

	ligandname = '../input/a01_resp'

	ligand = '../input/a01_resp.mol2'

	receptor = '../input/6td3_E_F.pdb'

	complex = '6td3_a01'

	neutral_expression = 'addions2 complex NA 2'

	tleap_input = make_tleap_input(ligandname, ligand, receptor, complex, neutral_expression)

	output = textwrap.dedent(f'''\
	source leaprc.protein.ff19SB
	source leaprc.phosaa19SB
	source leaprc.gaff
	loadamberparams ../input/a01_resp.frcmod
	loadamberparams frcmod.ionsjc_tip3p
	source leaprc.water.tip3p
	LIG = loadMol2 ../input/a01_resp.mol2
	receptor = loadPDB ../input/6td3_E_F.pdb
	complex = combine {{receptor LIG}}
	set default PBRadii mbondi2
	addions2 complex NA 2
	solvatebox complex TIP3PBOX 10.0
	savepdb complex 6td3_a01-box.pdb
	saveAmberParm complex 6td3_a01.parm7 6td3_a01.rst7
	quit
	''')	

	assert tleap_input == output


	os.chdir('../')


def test_get_neutral_expression():

	os.chdir('output')

	ligandname = '../input/a01_resp'

	ligand = '../input/a01_resp.mol2'

	receptor = '../input/6td3_E_F_protein.pdb'

	complex = '6td3_a01'

	neutral_expression = ''

	os.system('rm -rf tleap.in')
	tleap_input = make_tleap_input(ligandname, ligand, receptor, complex, neutral_expression)

	with open('tleap.in', 'w') as f:
		f.write(tleap_input)

	tleap_output = run(f'tleap -f tleap.in')

	print(tleap_output)

	expression = get_neutral_expression(tleap_output)

	print('EXPRESSION', expression)

	assert expression == 'addions2 complex CL 1'

	os.chdir('../')


def test_run_tleap():

	os.system('cp input/a01_resp.mol2 output/')

	os.system('cp input/6td3_E_F_protein.pdb output/')	

	os.chdir('output')

	ligand = 'a01_resp.mol2'

	receptor = '6td3_E_F_protein.pdb'

	complex = '6td3_a01'	

	tleap_input = run_tleap(ligand, receptor, complex)

	print(tleap_input)

	check = textwrap.dedent('''\
	source leaprc.protein.ff19SB
	source leaprc.phosaa19SB
	source leaprc.gaff
	loadamberparams a01_resp.frcmod
	loadamberparams frcmod.ionsjc_tip3p
	source leaprc.water.tip3p
	LIG = loadMol2 a01_resp.mol2
	receptor = loadPDB 6td3_E_F_protein.pdb
	complex = combine {receptor LIG}
	set default PBRadii mbondi2
	addions2 complex CL 1
	solvatebox complex TIP3PBOX 10.0
	savepdb complex 6td3_a01-box.pdb
	saveAmberParm complex 6td3_a01.parm7 6td3_a01.rst7
	quit
	''')

	print(check)

	assert tleap_input == check

	os.chdir('../')

