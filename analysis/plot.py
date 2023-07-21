#!/usr/bin/python3

import numpy as np
import pandas as pd
import pdb
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse

from sklearn.metrics import r2_score

parser = argparse.ArgumentParser(description='')

# Required arguments
parser.add_argument('-i','--input', help='Input File', required=True)
parser.add_argument('-o','--output', help='Output file', required=True)
parser.add_argument('-x','--xprop', help='x property to plot', required=True)
parser.add_argument('-y','--yprop', help='y property to plot', required=True)

# Optional arguments
parser.add_argument('-xt','--xtitle', help='Title for x axis', required=False)
parser.add_argument('-yt','--ytitle', help='Title for y axis', required=False)
parser.add_argument('-t','--title', help='Plot title', required=False)
parser.add_argument('-r2', action='store_true', help='Include coefficient of determination in plot title', required=False)
parser.add_argument('-s', '--select', help='Select data by column name and column value, e.g. -s "col_name 1", only works for int value', required=False)

args = parser.parse_args()

df = pd.read_csv(args.input, sep=';', index_col=0)

df.dropna(inplace=True, subset=[args.xprop, args.yprop])

fig, ax = plt.subplots()

if args.select is not None:

    col_name = args.select.split()[0]
    col_value = int(args.select.split()[1])

    df = df.loc[df[col_name] == col_value]

    print(f'Columns with {col_name} = {col_value} selected')


ax.scatter(df[args.xprop], df[args.yprop], s=1, color='red')

if args.xtitle:
    plt.xlabel(args.xtitle)
else:
    plt.xlabel(args.xprop)

if args.ytitle:
    plt.ylabel(args.ytitle)
else:
    plt.ylabel(args.yprop)

coefficient_of_determination = r2_score(df[args.xprop], df[args.yprop])

print(coefficient_of_determination)

if coefficient_of_determination < 0:
    coefficient_of_determination = 0

if args.r2:
    plt.title('R$\mathregular{^2}$ = %.2f' %coefficient_of_determination)

if args.title is not None:
    plt.title(args.title)

plt.savefig(args.output, dpi=600)

print(df)

