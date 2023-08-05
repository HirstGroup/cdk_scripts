import os
import sys

infile = open(sys.argv[1])

lig_list = []

for line in infile:
	lig1, lig2 = line.strip().split()
	print(lig1, lig2)

	#os.system(f'bash ~/cdk_scripts/mcl1/extractcomplex.sh {lig1} {lig2}')

	lig_list.append(lig1)

for i in lig_list:
	print(i, end=' ')
	