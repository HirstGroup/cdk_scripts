#!/usr/bin/python3

import argparse
import os
import sys
import textwrap

def calc_rmsd(parm, traj, output):
    """
    Calculate RMSD of MD

    Parameters
    ----------
    parm : str
        name of parm file
    traj : str
    	name of trajectory file
    output : str
		name of output file

    """

    string = textwrap.dedent(f'''\
    trajin {traj}
    rms ToFirst !@H= first out {output}
    ''')

    with open('rmsd.ptraj', 'w') as f:
    	f.write(string)

    print('line30', parm, traj, output)
    os.system(f'cpptraj {parm} rmsd.ptraj')


def calc_dist(parm, traj, output, mask1, mask2):
    """
    Calculate distance between mask1 and mask2 in MD

    Parameters
    ----------
    parm : str
        name of parm file
    traj : str
    	name of trajectory file
    output : str
		name of output file
	mask1 : str
		mask1 for distance
	mask2 : str
		mask2 for distance

    """

    string = textwrap.dedent(f'''\
    trajin {traj}
    distance dist {mask1} {mask2} out {output}
    ''')

    with open('dist.ptraj', 'w') as f:
    	f.write(string)

    print('line30', parm, traj, output)
    os.system(f'cpptraj {parm} dist.ptraj')


def get_atom_number_pdb(complex, LIG, atom_number_mol):
	"""
	Get atom number in full complex pdb from atom number in molecule

	Parameters
	----------
	complex : str
		name of complex
	LIG : str
		ligand molecule name
	atom_number_mol : str or int
		number of atom in ligand molecule

	Returns
	-------
	atom_number_pdb : str
		number of atom in complex pdb

	"""

	atom_number_mol = int(atom_number_mol)

	with open(f'{complex}-box.pdb') as f:

		n = 0

		for line in f:
			if LIG in line[17:20]:
				n += 1
				if n == atom_number_mol:
					atom_number_pdb = line.split()[1]
					return(atom_number_pdb)

	sys.exit(f'{complex} {LIG} {atom_number_mol} atom pdb not found')


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Calculate RMSD of MD')

	# Required arguments
	parser.add_argument('-c','--complex', help='Complex name', required=True)
	parser.add_argument('-m','--mask2', help='Mask2 to calculate distance', required=True)
	parser.add_argument('-r','--repeat', help='Repeat pattern of MD, e.g. _2, _3', default = '', required=False)

	# Optional arguments

	args = parser.parse_args()

	complex = args.complex

	parm = f'{complex}_strip.parm7'

	repeat = args.repeat

	traj = f'{complex}{repeat}_equi_cent_strip.nc'

	lig = complex.split('_')[1]

	LIG = lig.upper()

	mask1 = ':324@SG'

	atom_number_pdb = get_atom_number_pdb(complex, LIG, args.mask2)

	mask2 = f'@{atom_number_pdb}'

	calc_rmsd(parm, traj, f'{complex}{repeat}_equi_rmsd.dat')

	calc_dist(parm, traj, f'{complex}{repeat}_equi_dist.dat', mask1, mask2)



