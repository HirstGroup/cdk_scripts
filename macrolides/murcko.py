import argparse
import pandas as pd

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def get_basic_marcko_scaffold(smiles):
    return MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(smiles),includeChirality=False)


def get_bm_scaffold(smiles):
    try:
        scaffold = Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(mol=Chem.MolFromSmiles(smiles)))
    except Exception:
        print("Raise AtomValenceException, return basic Murcko Scaffold")
        scaffold = smiles
    return scaffold


def get_smiles_index(smiles_list):
    """
    Get index of unique smiles
    Used e.g. to get index of a Murcko scaffold
    """

    smiles_dict = dict()

    index_list = []

    for smi in smiles_list:
        if smi not in smiles_dict:
            smiles_dict[smi] = len(smiles_dict)
        index_list.append(smiles_dict[smi])

    return index_list


def main(input, output, smiles_col='SMILES'):

    df = pd.read_csv(input, sep=';')

    df['basic_murcko_scaffold'] = df[smiles_col].apply(get_basic_marcko_scaffold)
    df['BM_scaffold'] = df['basic_murcko_scaffold'].apply(get_bm_scaffold)
    df['BM_scaffold_index'] = get_smiles_index(df['BM_scaffold'])

    index_list = df['BM_scaffold_index'].to_list()

    df['BM_scaffold_index_count'] = [index_list.count(i) for i in index_list]

    df.to_csv(output, sep=';', index=False)

    return df


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calculate the Murcko scaffolds for a DataFrame')
    parser.add_argument('-i','--input', help='Input CSV File containing SMILES column', required=True)
    parser.add_argument('-o','--output', help='Output CSV File with calculated Murcko scaffolds', required=True)
    parser.add_argument('-c','--column', help='Name of column containing SMILES', default='SMILES', required=False)
    args = parser.parse_args()

    main(args.input, args.output, args.column)

