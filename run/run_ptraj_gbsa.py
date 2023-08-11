#!/usr/bin/env python

import argparse
import os
import sys

file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(f'{file_dir}/../')
from gbsa import antegbsa, gbsa
from ptraj import center_strip, delete_mds
from run import run

parser = argparse.ArgumentParser(description='Run a part of an MD (i.e. sequential MDs)')

# Required arguments
parser.add_argument('-c','--complex', help='Name of complexes to run', required=True)
parser.add_argument('-fp', '--first_part', type=int, help='First part to run MD, etc, e.g. 2, 3', required=False)
parser.add_argument('-lp', '--last_part', type=int, help='Last part to run MD, etc, e.g. 2, 3', required=True)

# Optional arguments
parser.add_argument('--ignore_errors', default='NO', help='Ignore errors when running system commands, i.e. continue and not crash, options YES/NO', required=False)
parser.add_argument('-r','--repeat', default='NO', help='Repeat pattern, e.g. _2, _3', required=False)
parser.add_argument('--test', default='NO', help='Do test run, i.e. very short MD', required=False)
parser.add_argument('-t', '--time', default='equi', help='Time pattern to run second MD, etc, e.g. equi, equi2', required=False)

args = parser.parse_args()

if args.last_part is None:
    args.last_part = args.first_part

if args.ignore_errors == 'YES':
    ignore_errors = True
else:
    ignore_errors = False

if args.repeat == 'NO':
    args.repeat = ''

if args.test == 'YES':
    print_cmd = True
else:
    print_cmd = False

for part in range(args.first_part, args.last_part + 1):

    center_strip(args.complex, ignore_errors=ignore_errors, print_cmd=print_cmd, repeat=args.repeat, time=f'{args.time}{part}')
    
    delete_mds(args.complex, repeat=args.repeat, time=f'{args.time}{part}')

    antegbsa(args.complex, part=part, print_cmd=print_cmd, repeat=args.repeat)
    
    gbsa(args.complex, part=part, print_cmd=print_cmd, repeat=args.repeat)