import filecmp

from murcko import *

def test_main():

	df = main('../tests/input/Nottingham_data15.csv', '../tests/output/murcko.csv', smiles_col='Smiles')

	print(df)

	assert list(df.head()['basic_murcko_scaffold']) == ['O=C(Nc1ccc(CN2CCNCC2)cc1)c1cccc(C#Cc2cnc3cccnn23)c1', 'O=C1Nc2ccccc2C1=Cc1ccccc1', 'O=C(Nc1cccc(Nc2nccc(-c3cccnc3)n2)c1)c1ccc(CN2CCNCC2)cc1', 'O=C(Nc1cccc(-n2ccnc2)c1)c1cccc(Nc2nccc(-c3cccnc3)n2)c1', 'c1ccc(Nc2nccc(-c3cnn4ncccc34)n2)cc1']
	
	assert list(df.head()['BM_scaffold']) == ['CC(CC1CCC(CC2CCCCC2)CC1)C1CCCC(CCC2CCC3CCCCC32)C1', 'CC1CC2CCCCC2C1CC1CCCCC1', 'CC(CC1CCCC(CC2CCCC(C3CCCCC3)C2)C1)C1CCC(CC2CCCCC2)CC1', 'CC(CC1CCCC(C2CCCC2)C1)C1CCCC(CC2CCCC(C3CCCCC3)C2)C1', 'C1CCC(CC2CCCC(C3CCC4CCCCC43)C2)CC1']


def test_get_smiles_index1():

	smiles_index = get_smiles_index(['C', 'C', 'CC', 'C', 'CCC', 'CC'])

	print(smiles_index)

	assert smiles_index == [0, 0, 1, 0, 2, 1]


def test_get_smiles_index2():

	df = pd.read_csv('../tests/input/murcko.csv', sep=';')

	smiles_list = df['BM_scaffold'].to_list()

	print(smiles_list)

	smiles_index = get_smiles_index(smiles_list)

	assert smiles_index == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 8, 6, 8, 8, 8, 8, 6, 10, 8, 11, 12, 11, 11, 13, 14, 8, 15, 6, 16, 16, 17, 17, 18, 16, 16, 16, 16, 18, 18, 16, 19, 20, 20, 20, 20, 18, 21, 22, 23, 24, 22, 24, 16, 20, 18, 22, 22, 25, 25, 18, 26, 18, 22, 22, 22, 27, 27, 16, 24, 24, 4, 4, 6, 28, 29, 6, 8, 30, 31, 32, 33, 32, 32, 8, 6, 34, 34, 35, 32, 6, 36, 34, 32, 32, 32, 32, 37, 38, 39, 40, 34, 32, 41, 42, 43, 6, 34, 39, 44, 6, 45, 46, 47, 6, 6, 48, 6, 39, 6, 34, 9, 49, 50, 50, 6, 6, 6, 6, 8, 8, 8, 8, 8, 51, 52, 53, 6, 6, 8, 52, 54, 52, 55, 56, 54, 8, 6, 55, 55, 8, 8, 34, 6, 57, 57, 12, 58, 59, 1, 58, 1, 8, 6, 60, 9, 61, 62, 6, 12, 6, 63, 63, 64, 6, 61, 6, 65, 66, 8, 8, 67, 8, 68, 6, 69, 70, 8, 12, 12, 6, 71, 8, 72, 15, 6, 73, 9, 6, 16, 16, 16, 16, 74, 75, 17, 18, 18, 76, 16, 17, 18, 77, 6, 22, 22, 22, 16, 18, 77, 78, 16, 21, 24, 18, 79, 80, 81, 81, 82, 27, 83, 24, 84, 85, 24, 86, 84, 21, 21, 18, 24, 84, 84, 24, 21, 21, 21, 18, 24, 24, 21, 21, 87, 21, 21, 21, 82, 88, 22, 22, 21, 63, 26, 24, 24, 24, 24, 26, 24, 89, 89, 90, 22, 91, 18, 84, 84, 24, 22, 6, 92, 18, 18, 22, 24, 24, 92, 93, 94, 18, 26, 26, 22, 95, 96, 79, 97, 98, 99, 100, 97, 101, 102, 103, 27, 93, 93, 104, 105, 98, 106, 107, 107, 108, 109, 18, 18, 110, 6, 109, 111, 112, 18, 109, 109, 18, 92, 24, 24, 56, 6, 24, 22, 22, 24, 16, 113, 89, 114, 109, 109, 16, 109, 22, 115, 20, 116, 16, 17, 117, 118, 119, 115, 113, 120, 120, 93, 120, 113, 113, 113, 113, 92, 22, 26, 118, 20, 121, 122, 24, 16, 24, 24, 6, 24, 123, 123, 124, 6, 26, 26, 16, 16, 16]

	assert filecmp.cmp('../tests/input/murcko.csv', '../tests/output/murcko.csv') is True

test_get_smiles_index2()
