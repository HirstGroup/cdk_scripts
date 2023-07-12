#!/usr/bin/env python

import argparse
import os
import subprocess
import sys
import textwrap


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


def antechamber(infile):
    """
    Run antechamber and parmchk2 on the input ligand structure

    Parameters
    ----------
    infile: ligand mol2 input file

    Returns
    -------
    None (generates output mol2 and frcmod files)
    """

    ligandname = os.path.splitext(infile)[0]    

    run(f'antechamber -i {infile} -fi mol2 -o {ligandname}_gas.mol2 -fo mol2 -c gas -dr no')

    run(f'parmchk2 -i {ligandname}_gas.mol2 -f mol2 -o {ligandname}_gas.frcmod')


def tleap(ligandname):
    """
    Run tleap

    Parameters
    ----------
    infile: input mol2 file

    Returns
    -------
    None (runs tleap)
    """

    tleap_input = textwrap.dedent(f'''\
    source leaprc.gaff
    loadamberparams {ligandname}.frcmod
    loadamberparams frcmod.ionsjc_tip3p
    complex = loadMol2 {ligandname}.mol2
    set default PBRadii mbondi2
    saveAmberParm complex {ligandname}.parm7 {ligandname}.rst7
    quit
    ''')

    with open('tleap.in', 'w') as f:
        f.write(tleap_input)

    run('tleap -f tleap.in')


def minimize(complex, bellymask):
    """
    Minimize ligand

    Parameters
    ----------
    complex: name of complex to look for input files (complex.parm7 and complex.rst7)

    Returns
    ------
    None (minimizes ligand and creates output files)
    """

    string = textwrap.dedent(f'''\
    Minimization
     &cntrl
       bellymask='{bellymask}'
       cut=999, 
       drms=0.001,
       DT=0.002, 
       ibelly=1,
       ig=-1,  
       igb=5,
       imin=1, 
       ioutfm=1, 
       irest=0, 
       maxcyc=5000, 
       ncyc=5000,
       NSTLIM=5000, 
       ntb=0, 
       ntmin=1, 
       ntp=0, 
       ntpr=100, 
       NTR=0, 
       ntwe=0,
       ntwr=5000, 
       ntwx=100, 
       ntx=1, 
     &end
    ''')

    with open('min.in', 'w') as f:
        f.write(string)

    run(f'sander -O -i min.in -c {complex}.rst7 -p {complex}.parm7 -o {complex}_min.out -r {complex}_min.rst7')


def make_bellymask_nonhydrogens(infile):
    """
    Generate bellymask for minimization of nonhydrogen atoms in input mol2 file

    Parameters
    ----------
    infile: input mol2 file

    Returns
    -------
    bellymask: string with atoms to minimize for Amber
    """

    atom_list = []

    with open(infile) as f:
        sel = False
        for line in f:
            if '@<TRIPOS>' in line: sel = False

            if sel:
                if line.split()[1][0] == 'H':
                    atom_list.append(line.split()[0])

            if '@<TRIPOS>ATOM' in line: sel = True

    bellymask = '@' + ','.join(atom_list)

    return bellymask


def convert_rst(parm, rst, output):
    """
    Convert rst7 to structure file, e.g. mol2 or pdb

    Parameters
    ----------
    parm: parm7 file
    rst: rst7 file
    output: output file

    Returns
    -------
    None (creates output file)
    """

    string = textwrap.dedent(f'''\
    trajin {rst}
    trajout {output} 
    ''')

    with open('convert_rst.ptraj', 'w') as f:
        f.write(string)

    run(f'cpptraj {parm} convert_rst.ptraj')


def optimize_hydrogens(infile, outfile):
    """
    Optimize hydrogens in mol2 structure keeping nonhydrogens fixed

    Parameters
    ----------
    infile: input mol2 file

    Returns
    -------
    None (generates output structure file)
    """

    ligandname = os.path.splitext(infile)[0]
    complex = ligandname + '_gas'

    antechamber(infile)

    tleap(complex)

    bellymask = make_bellymask_nonhydrogens(infile)

    #minimize(complex, bellymask)

    convert_rst(f'{complex}.parm7', f'{complex}_min.rst7', f'{complex}_min.pdb')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Minimize hydrogens in ligand mol2 file')
    parser.add_argument('-i','--input', help='Input mol2 file',required=True)
    parser.add_argument('-o','--output', help='Output structure file, e.g. mol2 or pdb',required=True)

    args = parser.parse_args()

    optimize_hydrogens(args.input, args.output)

