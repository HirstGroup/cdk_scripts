from generate_stereoisomers import *

def test1():

	smi = 'C=CC(=O)N1CCCC1C(=O)N1CCC(C1)c1cc2c([nH]1)nccc2c1cnn2c1cccn2'

	isomers = generate_stereoisomers(smi)

	output = ['C=CC(=O)N1CCC[C@@H]1C(=O)N1CC[C@@H](c2cc3c(-c4cnn5ncccc45)ccnc3[nH]2)C1', 'C=CC(=O)N1CCC[C@H]1C(=O)N1CC[C@@H](c2cc3c(-c4cnn5ncccc45)ccnc3[nH]2)C1', 'C=CC(=O)N1CCC[C@@H]1C(=O)N1CC[C@H](c2cc3c(-c4cnn5ncccc45)ccnc3[nH]2)C1', 'C=CC(=O)N1CCC[C@H]1C(=O)N1CC[C@H](c2cc3c(-c4cnn5ncccc45)ccnc3[nH]2)C1']

	assert isomers == output


def test2():

	smi = 'CN(C/C=C/C(=O)N1CC[C@@H]([C@@H](C1)C)C(=O)N1CCC(C1)c1cc2c([nH]1)nccc2c1cnn2c1cccn2)C'

	isomers = generate_stereoisomers(smi)

	output = ['C[C@@H]1CN(C(=O)/C=C/CN(C)C)CC[C@@H]1C(=O)N1CC[C@@H](c2cc3c(-c4cnn5ncccc45)ccnc3[nH]2)C1', 'C[C@@H]1CN(C(=O)/C=C/CN(C)C)CC[C@@H]1C(=O)N1CC[C@H](c2cc3c(-c4cnn5ncccc45)ccnc3[nH]2)C1']

	assert isomers == output


def test3():

	smi = 'CN(C/C=C/C(=O)N1CC[C@@H]([C@@H](C1)C)C(=O)N1CCC(C1)c1cc2c([nH]1)nccc2c1cnn2c1cccn2)C'

	isomers = generate_stereoisomers(smi)

	smiles_string = join_smiles(isomers)

	print(smiles_string)

	print(smiles_string.split('&'))

	assert isomers == smiles_string.split('&')

	
test3()
