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
parser.add_argument('-n','--n_repeats', type=int, help='Number of gbsa repeats', required=True)

# Optional arguments
parser.add_argument('-o','--output', help='Output CSV file name', required=False)
parser.add_argument('-p','--part', default='', help='Part pattern', required=False)

args = parser.parse_args()

part = args.part

df = pd.read_csv(args.input, sep=';')

df.dropna(inplace=True)

a = list(range(1,args.n_repeats+1))

df['gbsa_4_delta_total'] = df['gbsa_4_10_delta_total']

df[f'gbsa{part}_1_delta_total'] = df[f'gbsa{part}_delta_total']


def get_std(d_names, part=''):
	# Get standard deviation for dictionary of combinations vs iterations in GBSA data
	# d_names example {1: ['gbsa2_1', 'gbsa2_2', 'gbsa2_3'], 2: ['gbsa2_1_2', 'gbsa2_1_3', 'gbsa2_2_3'], 3: ['gbsa2_1_2_3']}

	for key, val in d_names.items():

		if len(val) == 1:
			df[f'gbsa{part}_{key}_std'] = 0
		else:			
			df[f'gbsa{part}_{key}_std'] = df[val].std(axis=1)


def get_avg(d, part=''):
	# Get average for dictionary d in GBSA data
	# Will create columns e.g. gbsa_1_2_3 for the average of gbsa_1, gbsa_2 and gbsa_3

	for key, val in d.items():

		df[key] = key

		# loop through combinations j will give tuples e.g. (1,2), (1,3)
		for j in val:

			k = [str(x) for x in j]

			name = '_'.join(k)

			print(name)

			df[f'gbsa{part}_{name}'] = 0

			# loop through tuples, e.g. loop through (1, 2) to get 1 and 2
			for i in j:
				df[f'gbsa{part}_{name}'] += df[f'gbsa{part}_{i}_delta_total']		

			df[f'gbsa{part}_{name}'] /= len(j)	


def make_dict_combinations(a):
	# Make dict mapping elements in list a to list of combinations (as tuples)
	# e.g. for 3 iterations d = {1: [(1,), (2,), (3,)], 2: [(1, 2), (1, 3), (2, 3)], 3: [(1, 2, 3)]}

	d = dict()

	for i in range(1, len(a)+1):

		l = []

		for j in list(combinations(a, i)):
			l.append(j)

		d[i] = l

	return d


def make_dict_names(d, part=''):
	# Make dictionary mapping iterations to list of names
	# e.g. for 3 iterations and part=2: {1: ['gbsa2_1', 'gbsa2_2', 'gbsa2_3'], 2: ['gbsa2_1_2', 'gbsa2_1_3', 'gbsa2_2_3'], 3: ['gbsa2_1_2_3']}

	d_names = dict()

	for key, val in d.items():

		l = []

		for j in val:

			k = [str(x) for x in j]

			name = f'gbsa{part}_' + '_'.join(k)

			l.append(name)

		d_names[key] = l

	return d_names


def plot_std(part=''):
	# Plot convergence graph of STD vs. No. simulations

	fig, ax = plt.subplots()

	for i in a:

		ax.scatter(df[i], df[f'gbsa{part}_{i}_std'])

	plt.xlabel('No. simulations')
	plt.ylabel('STD in energy (kcal/mol)')

	plt.xticks(a)

	plt.savefig(f'convergence{part}.pdf')

d = make_dict_combinations(a)

get_avg(d, part=part)

d_names = make_dict_names(d, part=part)
print(d_names)

get_std(d_names, part=part)

plot_std(part=part)

if args.output is not None:
	df.to_csv(args.output, sep=';', index=False)

print(df)

