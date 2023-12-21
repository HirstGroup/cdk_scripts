from rdkit import Chem

m = Chem.MolFromSmiles('CC')

m.GetSubstructMatch()