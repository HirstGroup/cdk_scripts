#!/usr/bin/env python

import argparse
import os 
import textwrap
import sys

file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(f'{file_dir}/../')
from run import run


def count_frames(parm, traj):
	"""
	Count number of frames in trajectory

	Parameters
	----------
	parm : str
		Name of parameter file
	traj : traj
		Name of trajectory file

	Returns
	-------
	n_frames : int
		Number of frames in trajectory
	"""

	string = textwrap.dedent(f'''\
	trajin {traj}
	run
	''')

	with open('count_frames.ptraj', 'w') as f:
		f.write(string)

	output = run(f'cpptraj {parm} count_frames.ptraj', verbose=False)

	n_frames = parse_count_frames_output(output)

	return n_frames


def parse_count_frames_output(count_frames_output):
	"""
	Parse count_frames output and return number of frames

	Parameters
	----------
	count_frames_output : str
		String containing cpptraj output from count_frames command

	Returns
	-------
	n_frames : int
		Number of frames
	"""

	for line in count_frames_output.splitlines():
		if 'frames and processed' in line:
			n_frames = int(line.split()[1])

			return n_frames


def delete_mds(complex, time='equi', repeat=''):
	"""
	Delete original MD in water and cent if cent_strip exists and has same number of frames

	Parameters
	----------
	complex : str
		Name of complex to look for {complex}.parm7 {complex}_{equi}.nc and {complex}_{equi}_cent.nc files

	Returns
	-------
	None (delete files)
	"""

	parm = f'{complex}.parm7'
	strip_parm = f'{complex}_strip.parm7'
	original_traj = f'{complex}{repeat}_{time}.nc'
	cent_traj = f'{complex}{repeat}_{time}_cent.nc'
	strip_traj = f'{complex}{repeat}_{time}_cent_strip.nc'

	if count_frames(parm, original_traj) == count_frames(parm, cent_traj):
		os.remove(original_traj)
		print(f'{original_traj} deleted')

	if count_frames(parm, cent_traj) == count_frames(strip_parm, strip_traj):
		os.remove(cent_traj)
		print(f'{cent_traj} deleted')


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Delete original MDs in water if cent_strip exists and has same number of frames')

	# Required arguments
	parser.add_argument('-i','--input', help='Input complex name', required=True)

	# Optional arguments
	parser.add_argument('-t','--time', default='equi', help='Time pattern to look for {complex}_time files', required=False)
	parser.add_argument('-r', '--repeat', default='', help='Repeat pattern, e.g. _2, _3', required=False)

	args = parser.parse_args()

	delete_mds(args.input, time=args.time, repeat=args.repeat)