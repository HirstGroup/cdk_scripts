import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem

input = sys.argv[1]

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


fps = [AllChem.GetMorganFingerprintAsBitVect(x,2,2048) for x in Chem.SmilesMolSupplier(input) if x is not None]

clusters=ClusterFps(fps,cutoff=0.4)

print(len(clusters))