import argparse
import os

parser = argparse.ArgumentParser(description='Run MD ptraj and gbsa')

# Required arguments
parser.add_argument('-c','--complex', help='Name of complex', required=True)

# Optional arguments
parser.add_argument('-p', '--part', default='', help='Part pattern to run second MD, etc, e.g. 2, 3', required=False)
parser.add_argument('-r','--repeat', default='', help='Repeat pattern, e.g. _2, _3', required=False)
parser.add_argument('-t','--test', default='NO', help='Run as a test', required=False)

args = parser.parse_args()

jobid=os.popen(f'sbatch --parsable ~/scripts/1gpu.sh bash ~/cdk_scripts/standard_md{args.part}.sh -c {args.complex} -r {args.repeat} --test {args.test}').read().splitlines()[0]

print(f'Submitted batch job {jobid}')

os.system(f'sbatch --dependency=afterok:{jobid} ~/scripts/1cpu.sh bash ~/cdk_scripts/run/run_ptraj_gbsa.sh -c {args.complex} -p {args.part} -r {args.repeat}')
