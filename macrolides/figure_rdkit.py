import argparse

parser = argparse.ArgumentParser(description='Make Figure of Cluster from RDKit output')
parser.add_argument('-i','--input', help='Input File Containing RDKit output (one tuple of cluster members per line)',required=True)
parser.add_argument('-s','--smiles', help='Aux File Containing smiles',required=True)
parser.add_argument('-o','--output', help='Output File',required=True)
parser.add_argument('-c', '--cluster', type=int, help='Cluster number to make Figure of',required=True)
args = parser.parse_args()

smifile = open(args.smiles)

smiles = []

for line in smifile:
	smiles.append(line.split()[0].strip())

infile = open(args.input)

clusters = []

for line in infile:
	clusters.append(eval(line))

outfile = open(args.output, 'w')

for i in clusters[args.cluster]:
	outfile.write(f'{smiles[i]}\n')
