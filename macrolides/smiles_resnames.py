import pandas as pd

def print_smiles_resnames_row(row):

	smi_list = row['stereoisomers_list'].split('&')
	resname_list = row['resname_list'].split('&')

	for smi, resname in zip(smi_list, resname_list):
		print(smi, resname)

df = pd.read_csv('Nottingham_data15.csv', sep=';')

print(df)

df.apply(print_smiles_resnames_row, axis=1)