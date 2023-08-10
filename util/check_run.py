#!/usr/bin/env python

import argparse
import os 
import textwrap
import sys

file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(f'{file_dir}/../')
from util.delete_md import count_frames
from run import run


def check_md(complex, frames, time, repeat, verbose=False):

	parm = f'{complex}.parm7'
	traj = f'{complex}{repeat}_{time}.nc'

	if not os.path.isfile(traj):
		if verbose: print(f'file {traj} not found')
		return 'NOK'

	if frames != count_frames(parm, traj, verbose=verbose):
		return 'NOK'

	return 'OK'

def check_md_cent(complex, frames, time, repeat, verbose=False):

	parm = f'{complex}.parm7'
	traj = f'{complex}{repeat}_{time}_cent.nc'

	if not os.path.isfile(traj):
		if verbose: print(f'file {traj} not found')
		return 'NOK'

	if frames != count_frames(parm, traj, verbose=verbose):
		return 'NOK'

	return 'OK'

def check_md_cent_strip(complex, frames, time, repeat, verbose=False):

	parm = f'{complex}_strip.parm7'
	traj = f'{complex}{repeat}_{time}_cent_strip.nc'

	if not os.path.isfile(traj):
		if verbose: print(f'file {traj} not found')
		return 'NOK'

	if frames != count_frames(parm, traj, verbose=False):
		return 'NOK'

	return 'OK'

def check_gbsa(complex, part, repeat, verbose=False):

	gbsa_file = f'gbsa{part}{repeat}/{complex}_gbsa{part}{repeat}.dat'

	if not os.path.isfile(gbsa_file):
		if verbose: print(f'file {gbsa_file} not found')
		return 'NOK'

	return 'OK'


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Check that MD ptraj and gbsa completed successfully')

	# Required arguments
	parser.add_argument('-i','--input', help='Input complex name', required=True)

	# Optional arguments
	parser.add_argument('-f','--frames', default=500, type=int, help='Number of frames in MD', required=False)
	parser.add_argument('-o','--output', help='Output file name containing results of check', required=False)
	parser.add_argument('-p','--part_list', nargs='+', default=[''], help='Parts to look for', required=False)
	parser.add_argument('-r', '--repeat_list', nargs='+', default=[''], help='Repeat pattern, e.g. _2, _3', required=False)
	parser.add_argument('-t','--time', nargs='+', default='equi', help='Time pattern to look for {complex}_time files', required=False)
	parser.add_argument('--verbose', action='store_true', help='Verbose output', required=False)

	args = parser.parse_args()

	complex = args.input
	frames = args.frames
	time = args.time
	verbose = args.verbose

	if args.output is not None:
		outfile = open(args.output, 'w')

		outfile.write(f'complex check_md check_md_cent check_md_cent_strip check_gbsa\n')

	for repeat in args.repeat_list:

		for part in args.part_list:

			if part == '1':
				part = ''

			check_md_out = check_md(complex, frames, time + part, repeat, verbose=verbose)

			check_md_cent_out = check_md_cent(complex, frames, time + part, repeat, verbose=verbose)			

			check_md_cent_strip_out = check_md_cent_strip(complex, frames, time + part, repeat, verbose=verbose)

			check_gbsa_out = check_gbsa(complex, part, repeat, verbose=verbose)

			print(complex, check_md_out, check_md_cent_out, check_md_cent_strip_out, check_gbsa_out)

			if args.output is not None:
				outfile.write(f'{complex} {check_md_out} {check_md_cent_out} {check_md_cent_strip_out} {check_gbsa_out}\n')




