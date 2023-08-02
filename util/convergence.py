#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from itertools import combinations

parser = argparse.ArgumentParser(description='Analyze convergence of GBSA')

# Required arguments
parser.add_argument('-i','--input', help='Input CSV file name', required=True)

args = parser.parse_args()

df = pd.read_csv(args.input, sep=';')

df.dropna(inplace=True)

a = [1,2,3,4]

df['gbsa_4_delta_total'] = df['gbsa_4_10_delta_total']

df['gbsa_1_delta_total'] = df['gbsa_delta_total']

print(df.columns)

def get_std(j):
	# Get standard deviation for tuple j in GBSA data
	pass

def get_avg(j):
	# Get average for tuple j in GBSA data

	k = [str(x) for x in j]

	name = '_'.join(k)

	print(name)

	df['gbsa_%s' %name] = 0

	for i in j:
		df['gbsa_%s' %name] += df['gbsa_%s_delta_total' %i]		

	df['gbsa_%s' %name] /= len(j)	

	return name

fig, ax = plt.subplots()

for i in range(1, len(a)+1):

	for j in list(combinations(a, i)):
		name = get_avg(j)

		df[name] = i

		ax.scatter(df[name], df['gbsa_%s' %name], s=1)

print(df)

plt.savefig('convergence.pdf')