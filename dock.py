#!/usr/bin/env python

import argparse
import numpy as np
import os
import pandas as pd
import sys

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

from gbsa import run

def generate_3d_sdf(smi, outfile_name):

    m = Chem.MolFromSmiles(smi)

    AllChem.Compute2DCoords(m)

    m = Chem.AddHs(m)

    AllChem.EmbedMolecule(m,randomSeed=0xf00d)

    print(Chem.MolToMolBlock(m),file=open(outfile_name,'w+'))


def protonate_mol2(infile_name, outfile_name):

    cmd = f'obabel -ph 7.4 -isdf {infile_name} -O {outfile_name}'

    print(cmd)

    os.system(cmd)


def protonate_mol2_row(row):

    for x, smi in enumerate(row['stereoisomers_list'].split('&')):

        ligname = '%s-%s' %(row['Row'], x+1)

        cmd = f'obabel -ph 7.4 -imol2 dock/6td3_{ligname}_dock_1.mol2 -O dock/6td3_{ligname}_dock_1_h.mol2'

        os.system(cmd)


def generate_mol2_row(row):

    for x, smi in enumerate(row['stereoisomers_list'].split('&')):

        generate_3d_sdf(smi, '%s-%s-3d.sdf' %(row['Row'], x+1))

        protonate_mol2('%s-%s-3d.sdf' %(row['Row'], x+1), '%s-%s.mol2' %(row['Row'], x+1))


def split_sdf(input, max_n=None):
    """
    Split SDF file into separate files, named as _n.sdf counting from 0

    Parameters
    ----------
    input: string
        name of input SDF file
    max_n: int
        maximum number of structures to extract


    Returns
    -------
    n_file: int
        total number of files created
    """

    name = input[:-4]

    with open(input) as infile:

        n = 0
        n_files = 0
        new = True

        for line in infile:
            if new:
                outfile = open('%s_%s.sdf' %(name, n_files), 'w')
                new = False
            outfile.write(line)
            if '$$$$' in line:
                n += 1
                n_files += 1
                new = True

                if n_files == max_n:
                    break

        outfile.close()

    return n_files

def dock(ligname):

    cmd = f'gnina -r 6td3_protein.pdb -l ligands/{ligname}.mol2 --autobox_ligand rc8.pdb --autobox_add 8 -o dockh/6td3_{ligname}_dock.sdf --log dockh/6td3_{ligname}_dock.out'

    os.system(cmd)


def dock_conf(ligname, folder_in, folder_out):
    """
    Dock molecule generating conformers for cycles using OpenEye (gnina does not generate conformers for cycles)

    Parameters
    ----------
    ligname: str
        name of ligand sdf file (without extension)
    folder_in: str
        folder where input file is located
    folder_out:
        folder where output files will be written to

    Returns
    -------
    None (writes docking sdf and out files for each conformer in output folder and then gets best conformation and writes it to *best.sdf and *best.out)
    """

    import generate_conformers_openeye

    # generate conformers with OpenEye
    generate_conformers_openeye.oe_conformer_generation2(f'{folder_in}/{ligname}', f'{folder_out}/{ligname}', tauto_sp23=False, torsion_drive=False, box_cen=None, save_mol2=True, save_conf_isomer_ids=True)

    # split SDF file
    n_files = split_sdf(f'{folder_out}/{ligname}_confs.sdf')

    # dock individual SDF files
    for i in range(n_files):
        run(f'gnina -r 6td3_protein.pdb -l {folder_out}/{ligname}_confs_{i}.sdf --autobox_ligand rc8.pdb --autobox_add 8 -o {folder_out}/{ligname}_confs_{i}_dock.sdf --log {folder_out}/{ligname}_confs_{i}_dock.out')


    # get lowest energy docking conformation
    score_list = []

    for i in range(n_files):
        score_list.append(parse_dock(f'{folder_out}/{ligname}_confs_{i}_dock.out'))

    index_min = np.argmin(score_list)

    os.system(f'cp {folder_out}/{ligname}_confs_{index_min}_dock.sdf {folder_out}/{ligname}_confs_dock_best.sdf')
    os.system(f'cp {folder_out}/{ligname}_confs_{index_min}_dock.out {folder_out}/{ligname}_confs_dock_best.out')

def write_gnina_output(score, outfile):
    """
    Write gnina style output file so that it can be parsed later

    Parameters
    ----------
    score: float
        docking score to write on output file
    outfile: str
        name of output file

    
    Returns
    -------
    None (writes output file)
    """

    with open(outfile, 'w') as f:

        f.write('-----+\n')
        f.write(f'1 {score}\n')


def dock_row(row):

    for ligname in row['resname_list'].split('&'):

        dock_conf(f'{ligname.lower()}', 'ligands', 'dock_conf')


def dock_parallel(df):

    with open('parallel.txt', 'w') as f:

        for index, row in df.iterrows():

            for x, smi in enumerate(row['stereoisomers_list'].split('&')):

                ligname = '%s-%s' %(row['Row'], x+1)

                f.write(f'gnina -r 6td3_protein.pdb -l {ligname}.mol2 --autobox_ligand rc8.pdb --autobox_add 8 -o 6td3_{ligname}_dock.mol2 --log 6td3_{ligname}_dock.out\n')

    os.system('parallel -j 4 -a parallel.txt')


def parse_dock(infile):
    """
    Parse docking score from gnina output

    Parameters
    ----------
    infile: str
        name of input file

    Returns
    -------
    score: float
        docking score

    """

    with open(infile) as f:

        sel = False

        for line in f:

            if sel:

                return float(line.split()[1])

            if '-----+' in line:

                sel = True    


def parse_dock_row(row):
    """
    Parse dock output row by row taking average of stereoisomers

    Parameters
    ----------
    row: pandas row

    Returns
    -------
    score_avg: float
        average of docking scores for all stereoisomers
    """

    if row['Covalent'] == True:
        return None

    #if pd.isnull(row['CDK12 Mean IC50 (uM)']):
    #    return None

    score_list = []

    for x, LIG in enumerate(row['resname_list'].split('&')):

        ligname = LIG.lower()

        score_list.append(parse_dock(f'dock_conf/{ligname}_confs_dock_best.out'))

    score_avg = np.mean(score_list)

    return score_avg


def parse_dock_row_ligname(row):
    """
    Parse dock output for a single enantiomer (single ligname)

    Parameters
    ----------
    row : pandas row

    Returns
    -------
    score : float
        Docking score
    """

    ligname = row['ligname']

    score = parse_dock(f'dock_conf/{ligname}_confs_dock_best.out')

    return score


def get_dock_conformation(row):
    """
    Extract first SDF structure from docking output
    into a new file for each row

    Parameters
    ----------
    row: pandas row

    Returns
    -------
    None (creates sdf file)
    """

    for x, LIG in enumerate(row['resname_list'].split('&')):

        ligname = LIG.lower()

        infile = 'dock_conf/%s_confs_dock_best.sdf' %ligname

        split_sdf(infile, max_n=1)


def get_first_mol2(infile, outfile):
    """
    Get first structure from mol2 file

    Parameters
    ----------
    infile: name of input mol2 file
    outfile: name of output mol2 file

    Returns
    -------
    None (creates output file)
    """

    lines = []

    with open(infile) as f:

        sel = False

        for line in f:
            lines.append(line)
            
            if sel is True:
                if '@<TRIPOS>MOLECULE' in line:
                    break


            if '@<TRIPOS>MOLECULE' in line:
                sel = True

    with open(outfile, 'w') as f:
        for line in lines:
            f.write(line)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run gbsa.py')

    # Required arguments
    parser.add_argument('-f','--function', help='Function name to run',required=True)
    parser.add_argument('-i','--input', help='Input file name, csv file separated by semicolon',required=True)

    # Optional argumenets
    parser.add_argument('-n','--name', help='Name of column for output file',required=False)
    parser.add_argument('-o','--output', help='Output file name',required=False)

    args = parser.parse_args()

    if args.input == args.output:
        os.system(f'cp {args.input} {args.input}.bk')

    df = pd.read_csv(args.input, sep=';')

    #df.sort_values(by='Row', inplace=True)

    #df = df.loc[df['Covalent'] == False]

    #df.dropna(inplace=True, subset=['CDK12 Mean IC50 (uM)'])

    #df.head(n=15).apply(generate_mol2_row, axis=1)

    if args.name is not None:
        name = args.name
    else:
        name = args.function

    df[name] = df.apply(eval(args.function), axis=1)

    # remove added column if function has None as output
    col = df[name].tolist()

    if len(set(col)) == 1 and col[0] is None:
        df.drop([name], axis=1, inplace=True)

    if args.output is not None:
        df.to_csv(args.output, sep=';', index=False)


        