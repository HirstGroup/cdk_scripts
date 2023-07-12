#!/usr/bin/env python

import numpy as np
import pandas as pd

df = pd.read_csv('Nottingham_data2.csv', sep=';')

df.drop_duplicates(inplace=True)

df.to_csv('Nottingham_data3.csv', sep=';', index=False)

df2 = df[df.duplicated(subset=['InChI'], keep=False)]

df2not = df[~df.duplicated(subset=['InChI'], keep=False)]

df2not.to_csv('Nottingham_data5.csv', sep=';', index=False)

#df2.sort_values(by='Original-ID', inplace=True)

df2.to_csv('nottingham_data_tidy_duplicated.csv', sep=';', index=False)

inchi_dict = dict()

for index, row in df2.iterrows():
    
    inchi = row['InChI']
    
    if inchi in inchi_dict:
        inchi_dict[inchi]['cdk12'].append(row['CDK12 Mean IC50 (uM)'])
        inchi_dict[inchi]['cdk9'].append(row['CDK9 Mean IC50 (uM)'])
        inchi_dict[inchi]['Original-ID'] += '-' + str(row['Original-ID'])
        inchi_dict[inchi]['ID'] += '-' + row['ID']
        inchi_dict[inchi]['Argenta Batch'] += '-' + row['Argenta Batch']
        if inchi_dict[inchi]['Smiles'] != row['Smiles']:
            print('SMILES', inchi)
    else:
        inchi_dict[inchi] = dict()
        inchi_dict[inchi]['cdk12'] = [row['CDK12 Mean IC50 (uM)']]
        inchi_dict[inchi]['cdk9'] = [row['CDK9 Mean IC50 (uM)']]
        inchi_dict[inchi]['ID'] = row['ID']
        inchi_dict[inchi]['Original-ID'] = str(row['Original-ID'])
        inchi_dict[inchi]['Argenta Batch'] = row['Argenta Batch']
        inchi_dict[inchi]['Smiles'] = row['Smiles']

inchi_col = []
cdk12_col = []
cdk12_std_col = []
cdk9_col =[]
cdk9_std_col = []
smiles_col = []
ID_col = []
original_ID_col = []
argenta_batch_col = []
smiles_col = []

for key, items in inchi_dict.items():
    inchi = key
    
    cdk12 = np.mean(items['cdk12'])
    cdk12_std = np.std(items['cdk12'])
    cdk9 = np.mean(items['cdk9'])
    cdk9_std = np.std(items['cdk9'])
    ID = items['ID']
    original_ID = items['Original-ID']
    argenta_batch = items['Argenta Batch']
    smiles = items['Smiles']
    
    inchi_col.append(inchi)
    cdk12_col.append(cdk12)
    cdk12_std_col.append(cdk12_std)
    cdk9_col.append(cdk9)
    cdk9_std_col.append(cdk9_std)
    ID_col.append(ID)
    original_ID_col.append(original_ID)
    argenta_batch_col.append(argenta_batch)
    smiles_col.append(smiles)
    
df3 = pd.DataFrame(list(zip(original_ID_col, ID_col, argenta_batch_col, cdk12_col, cdk9_col, smiles_col, inchi_col, cdk12_std_col, cdk9_std_col)), columns=['Original-ID', 'ID', 'Argenta Batch', 'CDK12 Mean IC50 (uM)', 'CDK9 Mean IC50 (uM)', 'Smiles', 'InChI', 'cdk12_std', 'cdk9_std'])

df3.to_csv('Nottingham_data4.csv', sep=';', index=False)

