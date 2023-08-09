#!/usr/bin/env python

import argparse
import os

parser = argparse.ArgumentParser(description='Run a part of an MD (i.e. sequential MDs)')

# Required arguments
parser.add_argument('-c','--complex', help='Name of complexes to run', required=True)

# Optional arguments
parser.add_argument('-p', '--part_list', type=int, nargs='+', help='List of parts to run MD, etc, e.g. 2, 3', required=False)
parser.add_argument('-r','--repeat', default='NO', help='Repeat pattern, e.g. _2, _3', required=False)
parser.add_argument('--test', default='NO', help='Do test run, i.e. very short MD', required=False)

args = parser.parse_args()

first_part = args.part_list[0]

cmd = f'sbatch --parsable ~/scripts/1gpu.sh python ~/cdk_scripts/standard_md_part.py -c {args.complex} -p {first_part} -r {args.repeat} --test {args.test}'

out = os.popen(cmd).read()

jobid = out.splitlines()[0]

print(f'Submitted batch job {jobid}')

os.system(f'sbatch --dependency=afterok:{jobid} ~/scripts/1cpu.sh python ~/cdk_scripts/run/run_ptraj_gbsa.sh -c {args.complex} -p {first_part} -r {args.repeat}')

if len(args.part_list) > 1:

    for part in args.part_list[1:]:

        out = os.popen(f'sbatch --dependency=afterok:{jobid} --parsable ~/scripts/1gpu.sh python ~/cdk_scripts/standard_md_part.py -c {args.complex} -p {part} -r {args.repeat} --test {args.test}').read()

        jobid = out.splitlines()[0]

        print(f'Submitted batch job {jobid}')

        os.system(f'sbatch --dependency=afterok:{jobid} ~/scripts/1cpu.sh python ~/cdk_scripts/run/run_ptraj_gbsa.sh -c {args.complex} -p {part} -r {args.repeat}')        