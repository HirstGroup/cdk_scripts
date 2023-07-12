#!/usr/bin/env python

import rdkit
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem import Draw

import pandas as pd

def generate_stereoisomers(smi):

	opts = StereoEnumerationOptions(unique=True)

	m = Chem.MolFromSmiles(smi)
	isomers = list(EnumerateStereoisomers(m, options=opts))
	
	smi_list = []

	for x in isomers:

		smi_list.append( Chem.MolToSmiles(x, isomericSmiles=True) )

	return smi_list


def join_smiles(smi_list):

	return '&'.join(smi_list)


def generate_stereoisomers_row(row):

	isomers = generate_stereoisomers(row['Smiles'])	

	isomers_string = join_smiles(isomers)

	return isomers_string

if __name__ == '__main__':

	df = pd.read_csv('Nottingham_data10.csv', sep=';')

	df.sort_values(by='Row', inplace=True)

	df['stereoisomers_list'] = df.apply(generate_stereoisomers_row, axis=1)

	df.to_csv('Nottingham_data11.csv', sep=';', index=False)
		