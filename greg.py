from rdkit import Chem

m = Chem.MolFromSmiles('C(=C(\F)/C)(/C)\Cl')

bnd = m.GetBondWithIdx(0)

bnd.GetStereo()

bnd.GetStereoAtoms()

print(list(bnd.GetStereoAtoms()))
