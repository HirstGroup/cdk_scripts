import numpy as np
import os
import pandas as pd
import sys

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

def generate_3d_sdf(smi, outfile_name):

	m = Chem.MolFromSmiles(smi)

	AllChem.Compute2DCoords(m)

	m = Chem.AddHs(m)

	AllChem.EmbedMolecule(m,randomSeed=0xf00d)

	print(Chem.MolToMolBlock(m),file=open(outfile_name,'w+'))


def protonate_mol2(infile_name, outfile_name):

	cmd = f'obabel -ph 7.4 -isdf {infile_name} -O {outfile_name}'

	print(cmd)

	os.system(cmd)


def generate_mol2_row(row):

	for x, smi in enumerate(row['stereoisomers_list'].split('&')):

		generate_3d_sdf(smi, '%s-%s-3d.sdf' %(row['Row'], x+1))

		protonate_mol2('%s-%s-3d.sdf' %(row['Row'], x+1), '%s-%s.mol2' %(row['Row'], x+1))


def dock(ligname):

	cmd = f'gnina -r 6td3_protein.pdb -l {ligname}.mol2 --autobox_ligand rc8.pdb --autobox_add 8 -o 6td3_{ligname}_dock.mol2 --log 6td3_{ligname}_dock100.out --exhaustiveness 100'

	os.system(cmd)


def dock_row(row):

	for x, smi in enumerate(row['stereoisomers_list'].split('&')):

		dock('%s-%s' %(row['Row'], x+1))


def dock_parallel(df):

	with open('parallel.txt', 'w') as f:

		for index, row in df.iterrows():

			for x, smi in enumerate(row['stereoisomers_list'].split('&')):

				ligname = '%s-%s' %(row['Row'], x+1)

				f.write(f'gnina -r 6td3_protein.pdb -l {ligname}.mol2 --autobox_ligand rc8.pdb --autobox_add 8 -o 6td3_{ligname}_dock.mol2 --log 6td3_{ligname}_dock.out\n')

	os.system('parallel -j 4 -a parallel.txt')


def parse_dock(row):

	if row['Covalent'] == True:
		return None

	if pd.isnull(row['CDK12 Mean IC50 (uM)']):
		return None

	print(row['CDK12 Mean IC50 (uM)'])

	score_list = []

	for x, smi in enumerate(row['stereoisomers_list'].split('&')):

		ligname = '%s-%s' %(row['Row'], x+1)

		with open(f'6td3_{ligname}_dock100.out') as f:

			sel = False

			for line in f:

				if sel:

					score_list.append( float(line.split()[1]) )

					break

				if '-----+' in line:

					sel = True


	return np.mean(score_list)


if __name__ == '__main__':

	df = pd.read_csv('Nottingham_data12.csv', sep=';')

	#df.sort_values(by='Row', inplace=True)

	#df = df.loc[df['Covalent'] == False]

	#df.dropna(inplace=True, subset=['CDK12 Mean IC50 (uM)'])

	#df.head(n=15).apply(generate_mol2_row, axis=1)

	#df.apply(dock_row, axis=1)

	df['Dock_score100'] = df.apply(parse_dock, axis=1)

	df.to_csv('Nottingham_data13.csv', sep=';', index=False)


		