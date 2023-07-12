#!/usr/bin/env python

import rdkit
from rdkit import Chem
from rdkit.Chem import Draw

import pandas as pd

df = pd.read_csv('Nottingham_data9.csv', sep=';')

print(df)

def remove_stereo(row):

    print(row['Smiles'])

    mol = Chem.MolFromSmiles(row['Smiles'])
    
    Chem.RemoveStereochemistry(mol) 

    return Chem.MolToSmiles(mol)
    
df['smi_flat'] = df.apply(remove_stereo, axis=1)

df.to_csv('Nottingham_data10.csv', sep=';', index=False)

