import argparse

parser = argparse.ArgumentParser(description='Make Figure of Cluster')
parser.add_argument('-i','--input', help='Input File',required=True)
parser.add_argument('-s','--smiles', help='Aux File',required=True)
parser.add_argument('-o','--output', help='Output File',required=True)
parser.add_argument('-c', '--cluster', type=int, help='Cluster number',required=True)
args = parser.parse_args()

smifile = open(args.smiles)

smiles = []

for line in smifile:
	smiles.append(line.strip())

infile = open(args.input)

n = 0
sel = False
ids = []

for line in infile:
	if sel:
		ids.extend([ int(x) for x in line.split()[1:] ])
		break

	if 'other members' in line:

		if n == args.cluster:
			sel = True
			ids = [int(line.split()[0])]

		n += 1

print(len(ids))

outfile = open(args.output, 'w')

for i in ids:
	outfile.write('%s\n' %(smiles[i]))

