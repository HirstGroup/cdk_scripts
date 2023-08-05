#!/usr/bin/env python

import argparse
import numpy as np
import sys

parser = argparse.ArgumentParser(description='Obtain unique resid from PDB')

# Required arguments
parser.add_argument('-i','--input', help='Input PDB file name',required=True)
parser.add_argument('-d','--dist', help='Distance to include in output file names',required=True)

args = parser.parse_args()

dist = args.dist

def ranges(p):
    q = sorted(p)
    i = 0
    for j in range(1,len(q)):
        if q[j] > 1+q[j-1]:
            yield q[i], q[j-1]
            i = j
    yield q[i], q[-1]


resid = set()

df = np.genfromtxt(args.input, names=None, skip_header=1, usecols=5, skip_footer=1, dtype=int)

print(df)

resid.update(set(df))

full = open('resid%sfull.txt' %dist, 'w')

for i in sorted( list(resid) ):
	full.write('%s ' %i)
full.write('\n')
	
resid_range = list(ranges(sorted(list(set(resid)))))

tex = open('resid%s.txt' %dist, 'w') 

for (x,y) in resid_range:
 print(x, y)
 if x!=y:
  tex.write('%s-%s,' %(x,y) )
 else:
  tex.write('%s,' %x)

tex.write('\n')
tex.close()

	
