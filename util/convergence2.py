#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from itertools import combinations

sys.exit('Work in progress to convert to functions')

def get_std(df, d_names, part=''):
	# Get standard deviation for dictionary of combinations vs iterations in GBSA data

	for key, val in d_names.items():

		if len(val) == 1:
			df[f'gbsa{part}_{key}_std'] = 0
		else:			
			df[f'gbsa{part}_{key}_std' %key] = df[val].std(axis=1)


def get_avg(d):
	# Get average for dictionary d in GBSA data

	sys.exit('get_avg not working')

	for key, val in d.items():

		df[key] = key

		for j in val:

			k = [str(x) for x in j]

			name = '_'.join(k)

			print(name)

			df['gbsa_%s' %name] = 0

			for i in j:
				df['gbsa_%s' %name] += df['gbsa_%s_delta_total' %i]		

			df['gbsa_%s' %name] /= len(j)	


def make_dict_combinations(a):
	# Make dict mapping elements in list a to list of combinations

	d = dict()

	for i in range(1, len(a)+1):

		l = []

		for j in list(combinations(a, i)):
			l.append(j)

		d[i] = l

	return d


def make_dict_names(d, part=''):
	# Make dictionary mapping iterations to list of names

	d_names = dict()

	for key, val in d.items():

		l = []

		for j in val:

			k = [str(x) for x in j]

			name = f'gbsa{part}_' + '_'.join(k)

			l.append(name)

		d_names[key] = l

	return d_names


def plot_std():
	# Plot convergence graph of STD vs. No. simulations

	fig, ax = plt.subplots()

	for i in a:

		ax.scatter(df[i], df['gbsa_%s_std' %i])

	plt.xlabel('No. simulations')
	plt.ylabel('STD in energy (kcal/mol)')

	plt.xticks(a)

	plt.savefig('convergence.pdf')


if __name__ == '__main__':

	arser = argparse.ArgumentParser(description='Analyze convergence of GBSA')

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

	a = list(range(1, args.n_repeats+1))

	df['gbsa_4_delta_total'] = df['gbsa_4_10_delta_total']

	df[f'gbsa{part}_1_delta_total'] = df[f'gbsa{part}_delta_total']

	d = make_dict_combinations(a)

	#get_avg(d)

	d_names = make_dict_names(d, part=args.part)
	print(d_names)

	get_std(df, d_names)

	plot_std()

	if args.output is not None:
		df.to_csv(args.output, sep=';', index=False)

	print(df)

