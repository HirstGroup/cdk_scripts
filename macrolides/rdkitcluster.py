import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem

import argparse

parser = argparse.ArgumentParser(description='Cluster molecules using RDKit Butina algorithm')
parser.add_argument('-i','--input', help='Input File Containing SMILES, one per line',required=True)
parser.add_argument('-o','--output', help='Output File of cluster tuples, one per line',required=True)
parser.add_argument('-c','--cutoff', help='Cutoff for clustering', type=float, required=True)
args = parser.parse_args()

#Define clustering setup
def ClusterFps(fps,cutoff=0.2):
    from rdkit import DataStructs
    from rdkit.ML.Cluster import Butina

    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
        dists.extend([1-x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True, reordering=True)
    return cs


fps = [AllChem.GetMorganFingerprintAsBitVect(x,2,2048) for x in Chem.SmilesMolSupplier(args.input) if x is not None]

clusters = ClusterFps(fps, cutoff=args.cutoff)

outfile = open(args.output, 'w')

for i in clusters:
    print(i)
    outfile.write(f'{i}\n')

print(len(clusters))