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
	None

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


def check_resp1_output(infile, row):

	with open(infile) as f:
	    lines = f.read().splitlines()
	    last_line = lines[-1]

	assert last_line.startswith(' Normal termination')