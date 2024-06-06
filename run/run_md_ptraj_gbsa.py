#!/usr/bin/env python

import argparse
import os

parser = argparse.ArgumentParser(description='Run MD ptraj and gbsa')

# Required arguments
parser.add_argument('-c','--complex', help='Name of complex', required=True)

# Optional arguments
parser.add_argument('-l','--ligandmask', help='Ligandmask', required=False)
parser.add_argument('-r','--repeat', default='NO', help='Repeat pattern, e.g. _2, _3', required=False)
parser.add_argument('-t','--test', default='NO', help='Run as a test', required=False)

args = parser.parse_args()

jobid=os.popen(f'sbatch --parsable ~/scripts/1gpu.sh bash ~/cdk_scripts/standard_md.sh -c {args.complex} -r {args.repeat} --test {args.test}').read().splitlines()[0]

print(f'Submitted batch job {jobid}')

if args.ligandmask is None:
	os.system(f'sbatch --dependency=afterok:{jobid} ~/scripts/1cpu.sh python ~/cdk_scripts/run/run_ptraj_gbsa.py -c {args.complex} -r {args.repeat}')
else:
	os.system(f'sbatch --dependency=afterok:{jobid} ~/scripts/1cpu.sh python ~/cdk_scripts/run/run_ptraj_gbsa.py -c {args.complex} -r {args.repeat} -l {args.ligandmask}')


