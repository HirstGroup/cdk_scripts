import pandas as pd
import sys,os

select = [int(x) for x in sys.argv[1].split()]

input = sys.argv[2]
output = sys.argv[3]

outfile = open(output, 'w')

df = pd.read_csv(input, sep=';')

df.sort_values(by='Row', inplace=True)

for index, row in df.iterrows():

    if row['Row'] in select:
        outfile.write('%s %s\n' %(row['Smiles'], row['Row']) )
        print('%s %s' %(row['Smiles'], row['Row']))
