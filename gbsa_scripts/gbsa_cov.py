#!/usr/bin/env python

import argparse
import os
import sys
import textwrap

file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(f'{file_dir}/../')
from run import run


def create_box(complex):
	"""
	Create box pdb for complex

	Parameters
	----------
	complex : str
		Complex name

	"""

	with open('box.ptraj', 'w') as f:
		f.write(f'trajin {complex}.rst7\n')
		f.write(f'trajout {complex}-box.pdb\n')

	run(f'cpptraj {complex}.parm7 box.ptraj')


def determine_strip_atoms(complex, cov_indeces, ligandmask):
	"""
	Determine atom numbers of ligand to strip

	Parameters
	----------
	complex : str
		Complex name
	cov_indeces : list of int
		List of indeces of atoms corresponding to covalent residue starting from zero
	ligandmask : str
		Mask of ligand

	Returns
	-------
	strip_atoms : list of int
		List of atom indeces starting from 1 to strip
	"""

	create_box(complex)

	ligand_atoms = []

	with open(f'{complex}-box.pdb') as f:

		for line in f:
			if line[17:20] == ligandmask:
				atom_number = int(line.split()[1])
				ligand_atoms.append(atom_number)

	first_atom = ligand_atoms[0]

	cov_atoms = [i + first_atom for i in cov_indeces]

	strip_atoms = [i for i in ligand_atoms if i not in cov_atoms]

	return strip_atoms


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run GBSA for covalent complex')

    # Required arguments
    parser.add_argument('-i','--input', help='Name of complexes to run functions on', required=True)

    # Optional arguments
    parser.add_argument('--cd', help='Name of directory to change into to run commands', required=False)
    parser.add_argument('--cov_indeces', nargs='+', type=int, help='List of indeces of atoms corresponding to covalent residue starting from zero', required=False)
    parser.add_argument('-f','--functions', nargs='+', help='Function names to run', required=False)
    parser.add_argument('-l','--ligandmask', help='Ligandmask', required=False)
    parser.add_argument('-p', '--part', default='', help='Part pattern to run second MD, etc, e.g. 2, 3', required=False)
    parser.add_argument('-r','--repeat', default='', help='Repeat pattern, e.g. _2, _3', required=False)
    parser.add_argument('--test', default='NO', help='Test run', required=False)

    args = parser.parse_args()

    complex = args.input

    determine_strip_atoms(complex, args.cov_indeces, args.ligandmask)