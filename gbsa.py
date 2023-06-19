import numpy as np
import os
import pandas as pd
import sys
import textwrap

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

def get_charge(infile):
	"""
	Get the total charge of a mol2 file

	Parameters
	----------
	infile: file name of mol2 file

	Returns
	-------
	charge: total charge of molecule

	"""

	infile = open(infile)

	charge = 0
	select = False

	for line in infile:
		if '@<TRIPOS>' in line: select = False
		if select: charge = charge + float(line.split()[-1])
		if '@<TRIPOS>ATOM' in line: select = True

	charge = int(np.around(charge))

	return charge


def create_resp1_file(infile, outfile, charge, cpu=1):
	"""
	Create resp1 file (optimization file in Gaussian)

	Parameters
	----------
	infile: name of molecule input file
	outfile: name of Gaussian output file
	charge: charge of molecule
	cpu: number of CPUs to use for Gaussian calculation

	Returns
	-------
	None (creates Gaussian file)

	"""

	string = textwrap.dedent(f'''\
	%nprocshared={cpu}
	#n HF 6-31G* opt

	{infile}

	{charge} 1
	''')

	with open(outfile, 'w') as f:
		f.write(string)

	os.system(f'obabel -imol2 {infile} -oxyz | tail -n+3 >> {outfile}')

	with open(outfile, 'a') as f:
		f.write('\n')


def check_resp1_output(infile, inchikey):
	"""
	Check output of resp1 calculation

	Parameters
	----------
	infile: name of Gaussian output file
	inchikey: inchikey of starting structure

	Returns
	-------
	check: 'OK' or 'NOK' string
	"""

	with open(infile) as f:
	    lines = f.read().splitlines()
	    last_line = lines[-1]

	output = os.popen(f'obabel -ig09 {infile} -omol2 | obabel -imol2 -oinchikey').read().split('\n')[0]

	check = 'OK'

	if inchikey != output:
		check = 'NOK'

	if 'Normal termination' not in last_line:
		check = 'NOK'

	return check

def get_inchikey(infile):
	"""
	Get inchikey for structure using obabel

	Parameters
	----------
	infile: structure input file

	Returns
	-------
	inchikey: inchikey of structure
	"""

	inchikey = os.popen(f'obabel {infile} -oinchikey').read().split('\n')[0]	

	return inchikey


def create_resp2_file(infile, outfile, charge, cpu=1):
	"""
	Create resp2 file (ESP electrostatic potential file in Gaussian)

	Parameters
	----------
	infile: name of Gaussian resp1 output file
	outfile: name of Gaussian output file
	charge: charge of molecule
	cpu: number of CPUs to use for Gaussian calculation

	Returns
	-------
	None (creates Gaussian file)

	"""

	string = textwrap.dedent(f'''\
	%nprocshared={cpu}
	#n HF 6-31G* POP(MK) IOP(6/33=2) scf=direct

	{infile}

	{charge} 1
	''')

	with open(outfile, 'w') as f:
		f.write(string)

	os.system(f'obabel -ig09 {infile} -oxyz | tail -n+3 >> {outfile}')

	with open(outfile, 'a') as f:
		f.write('\n')


def create_resp3_file(infile, outfile1, outfile2, auxfile, resname):
	"""
	Create mol2 file with resp charges with antechamber from Gaussian ESP file with original coordinates

	Parameters
	----------
	infile: output ESP file from Gaussian
	outfile: mol2 output file
	auxfile: mol2 file with original coordinates

	Returns
	-------
	None (creates mol2 file)
	"""

	with open(infile) as f:
	    lines = f.read().splitlines()
	    last_line = lines[-1]

	assert 'Normal termination' in last_line

	os.system(f'antechamber -i {infile} -fi gout -gv 1 -o {outfile1} -fo mol2 -c resp -rn {resname}')

	os.system(f'antechamber -i {outfile1} -fi mol2 -o {outfile2} -fo mol2 -a {auxfile} -ao crd')

	os.system('rm ANTECHAMBER* ATOMTYPE.INF esout punch qout QOUT')


def check_resp3_file(infile, outfile, inchikey, charge):
	"""
	Check that resp mol2 file is correct, in terms of charge and inchikey

	Parameters
	----------
	infile: input mol2 resp file
	inchikey: inchikey of structure
	charge: charge of structure

	Returns
	-------
	check: string 'OK' or 'NOK'
	"""

	sys.exit('Not working because obabel gives different can, inchi and inchikey for structures')

	check = 'OK'

	os.system(f'antechamber -imol2 {infile} -fi mol2 -o ')

	output = os.popen(f'obabel {outfile} -oinchikey').read().split('\n')[0]

	if inchikey != output:
		check = 'NOK'

	if get_charge(infile) != charge:
		check = 'NOK'

	return check




	
