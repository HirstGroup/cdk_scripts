#!/usr/bin/env python

import argparse
import os
import sys

file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(f'{file_dir}/../')

from parse_gbsa_decomp import *

parser = argparse.ArgumentParser(description='Parse GBSA decomp results')

# Required arguments
parser.add_argument('-c','--complex', help='Name of complex to run', required=True)

args = parser.parse_args()

complex = args.complex

df = parse_gbsa_decomp(f'{complex}/gbsa_decomp/{complex}_gbsa_decomp.dat')

input_list = [f'{complex}/gbsa_decomp/{complex}_gbsa_decomp.dat']

for i in range(2,6):
	input_list.append( f'{complex}/gbsa_decomp_{i}/{complex}_gbsa_decomp_{i}.dat' )

df = average_decomp(input_list, output=f'{complex}/{complex}_gbsa_decomp_avg.csv')

print(df)
