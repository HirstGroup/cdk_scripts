#!/usr/bin/env python

import os
import sys

from rdkit import Chem

import argparse

import pandas as pd
 

def pattern(row):

	x = row['Smiles']
	idx = row['Row']

	m = Chem.MolFromSmiles(x)
	match = False
	for pattern_i in pattern_list_objects:
		if m.HasSubstructMatch(pattern_i):
			print(idx, x)
			match = True

	return match


def has_pattern(infile, pattern):

	m = Chem.MolFromMolFile(infile)

	pattern_m = Chem.MolFromSmarts(pattern)

	return m.HasSubstructMatch(pattern_m)


def get_substruct(infile, pattern):
	
	pattern_m = Chem.MolFromSmarts(pattern)

	m = Chem.MolFromMolFile(infile)
	
	atom_list = m.GetSubstructMatch(pattern_m)

	return atom_list


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='This script splits SDF file')
	parser.add_argument('-i','--input', help='Input file name', required=True)
	parser.add_argument('-p','--pattern_list', nargs='+', help='Smarts pattern list space separated, use individual quotation marks for each pattern', required=False)
	parser.add_argument('-n','--pattern_name', help='Pattern Name', required=False)
	parser.add_argument('-o','--output', help='Output file name', required=True)
	args = parser.parse_args()

	if args.pattern_list is not None:

		pattern_list = args.pattern_list

		pattern_name = args.pattern_name

	else:

		# covalent pattern list
		pattern_list = ['[#6]=[#6]-[#6]=O', 'Cl[#6]-[#6]=O', 'O=[#6]C#C']

		pattern_name = 'covalent'

#df[pattern_name] = df.apply(pattern, axis=1)

#df.to_csv(args.output, sep=';', index=False)
