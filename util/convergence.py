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

a = list(range(1,11))

df['gbsa_4_delta_total'] = df['gbsa_4_10_delta_total']

df['gbsa_1_delta_total'] = df['gbsa_delta_total']


for i in range(5,11):
	df['gbsa_%s_delta_total' %i] = df['gbsa_1_delta_total']

def get_std(d_names):
	# Get standard deviation for dictionary of combinations vs iterations in GBSA data

	for key, val in d_names.items():

		if len(val) == 1:
			df['gbsa_%s_std' %key] = 0
		else:			
			df['gbsa_%s_std' %key] = df[val].std(axis=1)


def get_avg(d):
	# Get average for dictionary d in GBSA data

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


def make_dict_names(d):
	# Make dictionary mapping iterations to list of names

	d_names = dict()

	for key, val in d.items():

		l = []

		for j in val:

			k = [str(x) for x in j]

			name = 'gbsa_' + '_'.join(k)

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

d = make_dict_combinations(a)

get_avg(d)

d_names = make_dict_names(d)
print(d_names)

get_std(d_names)

plot_std()

print(df)

