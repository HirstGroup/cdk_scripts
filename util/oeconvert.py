import os
import sys

from openeye import oechem

ifs = oechem.oemolistream()
ifs.open(sys.argv[1])

ofs = oechem.oemolostream()
ofs.open(sys.argv[2])

ofs.SetFormat(oechem.OEFormat_MOL2)

for mol in ifs.GetOEGraphMols():
    oechem.OEWriteMolecule(ofs, mol)