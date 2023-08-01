#!/usr/bin/env python

import argparse 

parser = argparse.ArgumentParser(description='Find amino acid sequence on PDB file and return resid')
parser.add_argument('-i','--input', help='Input PDB File',required=True)
parser.add_argument('-g','--grep', nargs='+', help='Space separated list of contiguous amino acids to grep',required=True)

args = parser.parse_args()

grep = args.grep

f = open(args.input)

sequence = dict()

for line in f:
	if not line.startswith('ATOM'):
		continue

	resname = line[17:20]
	resid = int(line[22:27].split()[0])

	sequence[resid] = resname

all_resid_list = []
resid_list = []

for key, val in sequence.items():

	x = len(resid_list)
	if x == len(grep):
		all_resid_list.append(resid_list)
		resid_list = []
		x = 0

	if val == grep[x]:
		resid_list.append(key)
	else:
		resid_list = []

print(args.input)

for i in all_resid_list:
	for j in i:
		print(j, end=' ')
	print()
