#!/usr/bin/env python

import argparse
import os
import sys

file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(f'{file_dir}/../')
from run import run

parser = argparse.ArgumentParser(description='Run a part of an MD (i.e. sequential MDs)')

# Required arguments
parser.add_argument('-c','--complex', help='Name of complexes to run', required=True)

# Optional arguments
parser.add_argument('-fp', '--first_part', type=int, help='First part to run MD, etc, e.g. 2, 3', required=False)
parser.add_argument('-lp', '--last_part', type=int, help='Last part to run MD, etc, e.g. 2, 3', required=False)
parser.add_argument('-r','--repeat', default='NO', help='Repeat pattern, e.g. _2, _3', required=False)
parser.add_argument('--test', default='NO', help='Do test run, i.e. very short MD', required=False)

args = parser.parse_args()

first_part = args.first_part

if args.test == 'YES':
    print_cmd = True
else:
    print_cmd = False

cmd = f'sbatch --parsable ~/scripts/1gpu.sh ~/cdk_scripts/standard_md_part.py -c {args.complex} -p {first_part} -r {args.repeat} --test {args.test}'

out = run(cmd, print_cmd=print_cmd)

jobid = out.splitlines()[0]

print(f'Submitted batch job {jobid}')

run(f'sbatch --dependency=afterok:{jobid} ~/scripts/1cpu.sh ~/cdk_scripts/run/run_ptraj_gbsa.sh -c {args.complex} -p {first_part} -r {args.repeat}', print_cmd=print_cmd)

if args.last_part is not None:

    for part in range(args.first_part + 1, args.last_part + 1):

        cmd = f'sbatch --dependency=afterok:{jobid} --parsable ~/scripts/1gpu.sh ~/cdk_scripts/standard_md_part.py -c {args.complex} -p {part} -r {args.repeat} --test {args.test}'

        out = run(cmd, print_cmd=print_cmd)

        jobid = out.splitlines()[0]

        print(f'Submitted batch job {jobid}')

        run(f'sbatch --dependency=afterok:{jobid} ~/scripts/1cpu.sh ~/cdk_scripts/run/run_ptraj_gbsa.sh -c {args.complex} -p {part} -r {args.repeat}', print_cmd=print_cmd)        