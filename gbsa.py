import numpy as np
import os
import pandas as pd
import subprocess
import sys
import textwrap

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

def run(cmd):
    """
    Run a command using subprocess and make code crash if there's an error
    and print error on console

    Parameters
    ----------
    cmd: command to run

    Returns
    -------
    None (runs command)
    """

    subprocess.run(cmd, shell=True, check=True, stderr=subprocess.STDOUT)


def get_charge(infile):
    """
    Get the total charge of a mol2 file

    Parameters
    ----------
    infile: file name of mol2 file

    Returns
    -------
    charge: total charge of molecule

    """

    infile = open(infile)

    charge_list = []
    select = False

    for line in infile:
        if '@<TRIPOS>' in line: select = False
        if select: 
            charge_list.append(float(line.split()[-1]))
        if '@<TRIPOS>ATOM' in line: select = True

    if len(charge_list) == 0:
        sys.exit('No charges read')

    charge = int(np.around(sum(charge_list)))

    return charge


def create_resp1_file(infile, outfile, charge, cpu=1):
    """
    Create resp1 file (optimization file in Gaussian)

    Parameters
    ----------
    infile: name of molecule input file
    outfile: name of Gaussian output file
    charge: charge of molecule
    cpu: number of CPUs to use for Gaussian calculation

    Returns
    -------
    None (creates Gaussian file)

    """

    string = textwrap.dedent(f'''\
    %nprocshared={cpu}
    #n HF 6-31G* opt

    {infile}

    {charge} 1
    ''')

    with open(outfile, 'w') as f:
        f.write(string)

    os.system(f'obabel {infile} -oxyz | tail -n+3 >> {outfile}')

    with open(outfile, 'a') as f:
        f.write('\n')


def check_resp1_output(infile, inchikey):
    """
    Check output of resp1 calculation

    Parameters
    ----------
    infile: name of Gaussian output file
    inchikey: inchikey of starting structure

    Returns
    -------
    check: 'OK' or 'NOK' string
    """

    with open(infile) as f:
        lines = f.read().splitlines()
        last_line = lines[-1]

    output = os.popen(f'obabel -ig09 {infile} -omol2 | obabel -imol2 -oinchikey').read().split('\n')[0]

    check = 'OK'

    if inchikey != output:
        check = 'NOK'

    if 'Normal termination' not in last_line:
        check = 'NOK'

    return check

def get_inchikey(infile):
    """
    Get inchikey for structure using obabel

    Parameters
    ----------
    infile: structure input file

    Returns
    -------
    inchikey: inchikey of structure
    """

    inchikey = os.popen(f'obabel {infile} -oinchikey').read().split('\n')[0]    

    return inchikey


def create_resp2_file(infile, outfile, charge, cpu=1):
    """
    Create resp2 file (ESP electrostatic potential file in Gaussian)

    Parameters
    ----------
    infile: name of Gaussian resp1 output file
    outfile: name of Gaussian output file
    charge: charge of molecule
    cpu: number of CPUs to use for Gaussian calculation

    Returns
    -------
    None (creates Gaussian file)

    """

    string = textwrap.dedent(f'''\
    %nprocshared={cpu}
    #n HF 6-31G* POP(MK) IOP(6/33=2) scf=direct

    {infile}

    {charge} 1
    ''')

    with open(outfile, 'w') as f:
        f.write(string)

    os.system(f'obabel -ig09 {infile} -oxyz | tail -n+3 >> {outfile}')

    with open(outfile, 'a') as f:
        f.write('\n')


def create_resp3_file(infile, outfile1, outfile2, auxfile, resname):
    """
    Create mol2 file with resp charges with antechamber from Gaussian ESP file with original coordinates

    Parameters
    ----------
    infile: output ESP file from Gaussian
    outfile: mol2 output file
    auxfile: mol2 file with original coordinates

    Returns
    -------
    None (creates mol2 file)
    """

    with open(infile) as f:
        lines = f.read().splitlines()
        last_line = lines[-1]

    assert 'Normal termination' in last_line

    run(f'antechamber -i {infile} -fi gout -gv 1 -o {outfile1} -fo mol2 -c resp -rn {resname} -dr no')

    run(f'antechamber -i {outfile1} -fi mol2 -o {outfile2} -fo mol2 -a {auxfile} -fa mol2 -ao crd -dr no')

    os.system('rm ANTECHAMBER* ATOMTYPE.INF esout punch qout QOUT')


def check_resp3_file(infile, outfile, inchikey, charge):
    """
    Check that resp mol2 file is correct, in terms of charge and inchikey

    Parameters
    ----------
    infile: input mol2 resp file
    inchikey: inchikey of structure
    charge: charge of structure

    Returns
    -------
    check: string 'OK' or 'NOK'
    """

    sys.exit('Not working because obabel gives different can, inchi and inchikey for structures')

    check = 'OK'

    os.system(f'antechamber -imol2 {infile} -fi mol2 -o ')

    output = os.popen(f'obabel {outfile} -oinchikey').read().split('\n')[0]

    if inchikey != output:
        check = 'NOK'

    if get_charge(infile) != charge:
        check = 'NOK'

    return check


def run_tleap(ligand, receptor, complex):
    """
    Run tleap to make protein ligand complex

    Parameters
    ----------
    ligand: ligand file in mol2 format
    receptor: receptor file in pdb format
    complex: name of complex

    Returns
    -------
    None (complex parm7 and rst7 files created)
    """

    ligandname = os.path.splitext(ligand)[0]

    run(f'parmchk2 -i {ligand} -f mol2 -o {ligandname}.frcmod')

    tleap_input = make_tleap_input(ligandname, ligand, receptor, complex, '')
    
    with open('tleap.in', 'w') as outfile:
        outfile.write(tleap_input)

    tleap_output = os.popen(f'tleap -f tleap.in').read()

    neutral_expression = get_neutral_expression(tleap_output)

    if neutral_expression != '':
        tleap_input = make_tleap_input(ligandname, ligand, receptor, complex, neutral_expression)
    
        with open('tleap.in', 'w') as outfile:
            outfile.write(tleap_input)

        run(f'tleap -f tleap.in') 


def make_tleap_input(ligandname, ligand, receptor, complex, neutral_expression):
    """
    Make tleap input file

    Parameters
    ----------
    ligandname: root of ligand file name
    ligand: ligand file in mol2 format
    receptor: receptor file in pdb format
    complex: name of complex
    neutral_expression: expression used to neutralize complex
    outfolder: name of folder to write output to

    Returns
    -------
    tleap_input: string with tleap input file
    """

    tleap_input = textwrap.dedent(f'''\
    source leaprc.protein.ff19SB
    source leaprc.phosaa19SB
    source leaprc.gaff
    loadamberparams {ligandname}.frcmod
    loadamberparams frcmod.ionsjc_tip3p
    source leaprc.water.tip3p
    LIG = loadMol2 {ligand}
    receptor = loadPDB {receptor}
    complex = combine {{receptor LIG}}
    set default PBRadii mbondi2
    {neutral_expression}
    solvatebox complex TIP3PBOX 10.0
    savepdb complex {complex}-box.pdb
    saveAmberParm complex {complex}.parm7 {complex}.rst7
    quit
    ''')

    return tleap_input



def get_neutral_expression(tleap_output):
    """
    Get charge of complex and hence expression to neutralize complex in tleap

    Parameters
    ----------
    tleap_output: tleap_output string that will be parsed

    Returns
    -------
    expression: command to neutralize complex in tleap
    """

    charge = 0

    for line in tleap_output:
        if 'WARNING: The unperturbed charge of the unit:' in line: charge = float(line.split()[7])

    charge = int(np.round(charge))

    if charge < 0: 
        return 'addions2 complex NA %s' %-charge
    elif charge > 0: 
        return 'addions2 complex CL %s' %charge
    else:
        return ''

