import rdkit
from rdkit import Chem
from rdkit.Chem import Draw

import pandas as pd

df = pd.read_csv('Nottingham_data6.csv', sep=';')

print(df)

def chiral(row):

    print(row['Smiles'])

    mol = Chem.MolFromSmiles(row['Smiles'])
    
    chiral_list = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    
    return len(chiral_list)

df['n_chiral'] = df.apply(chiral, axis=1)

def chiral_assigned(row):

    print(row['Smiles'])

    mol = Chem.MolFromSmiles(row['Smiles'])
    
    chiral_list = Chem.FindMolChiralCenters(mol, includeUnassigned=False, useLegacyImplementation=False)
    
    return len(chiral_list)

df['n_chiral_assigned'] = df.apply(chiral_assigned, axis=1)

def double_bond(row):

    mol = Chem.MolFromSmiles(row['Smiles'])

    print(rdkit.Chem.rdmolops.FindPotentialStereoBonds(mol))

df['double_bond'] = df.apply(double_bond, axis=1)

df.to_csv('Nottingham_data7.csv', sep=';', index=False)

