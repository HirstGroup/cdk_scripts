#!/usr/bin/env python

import os
import pandas as pd
import sys

df = pd.read_csv('../Nottingham_data10.csv', sep=';')

for index, row in df.iterrows():

	smiles = row['Smiles']

	id = row['Row']

	#with open(f'{id}.smi', 'w') as f:
#		f.write(f'{smiles} {id}')

	os.system(f'obabel -:"{smiles}" -O {id}.svg')