#!/usr/bin/env python

import argparse
import numpy as np
import os
import sys
import textwrap

file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(f'{file_dir}/../')
from util.delete_md import count_frames
from run import run


def get_amber_range(resid_list, range_file):
    """
    Convert resid_list to amber format
    """

    range_tuple = list(get_range_tuple(resid_list))

    print(range_tuple)

    amber_range = ''

    with open(range_file, 'w') as f:
        for (x,y) in range_tuple:
            if x!=y:
                amber_range += f'{x}-{y},'
                f.write(f'{x} to {y} ')
            else:
                amber_range += f'{x},'
                f.write(f'{x} ')
        f.write('\n')

    return amber_range


def get_range_tuple(resid_list):
    """
    Get list of tuples describing range
    where tuple (a, b) means range goes from a to b inclusive
    """
    q = sorted(list(resid_list))
    i = 0
    for j in range(1,len(q)):
        if q[j] > 1+q[j-1]:
            yield q[i], q[j-1]
            i = j
    yield q[i], q[-1]


def get_bs_mask(frames, mask, parm, traj, exclude=None, interval=1):
    """
    Get mask of residues in binding site
    """

    if not os.path.exists('mask'):
        os.makedirs('mask')

    os.chdir('mask')

    string = textwrap.dedent(f'''\
    trajin ../{traj} 1 last {interval}
    mask {mask}<:4.0 maskpdb resid_4a.pdb
    ''')

    with open('mask.ptraj', 'w') as f:
        f.write(string)

    run(f'cpptraj ../{parm} mask.ptraj', print_cmd=True)

    resid_set = set()

    for frame in range(1,frames+1):
        df = np.genfromtxt(f'resid_4a.pdb.{frames}', names=None, skip_header=1, usecols=4, skip_footer=1, dtype=int)
        resid_set.update(set(df))

    if exclude is not None:
        resid_set = [i for i in resid_set if i not in exclude]
    
    amber_range = get_amber_range(resid_set, 'resid_vmd.txt')

    os.chdir('../')

    return amber_range


def cluster_bs(complex, epsilon=1.0, exclude=None, folder_name='', interval=1, mask=None, method='dbscan', parm=None, repeat='', time='equi', traj=None):
    """
    Cluster trajectory by distance to ligand in binding site
    """

    if mask is None:
        mask = ':' + complex.upper()

    if parm is None:
        parm = f'../{complex}_strip.parm7'
    else:
        parm = f'../{parm}'
    if traj is None:
        traj = f'../{complex}{repeat}_{time}_cent_strip.nc'
    else:
        traj = f'../{traj}'

    folder = f'cluster{folder_name}{repeat}_{time}_{interval}_{epsilon}'

    os.system(f'rm -rf {folder}')

    os.makedirs(folder)

    os.chdir(folder)

    counted_frames = count_frames(parm, traj, verbose=True) 

    frames = counted_frames // interval

    if counted_frames % interval > 0:
        frames += 1

    bs_mask = get_bs_mask(frames, mask, parm, traj, exclude=exclude, interval=interval)

    if method == 'dbscan':
        string = textwrap.dedent(f'''\
        trajin {traj} 1 last {interval}
        cluster {complex} \
        dbscan minpoints 2 epsilon {epsilon} \
        rms {mask},:{bs_mask} nofit \
        sieve 1 random \
        pairdist pairdist.dat \
        out cnumvtime.dat \
        sil Sil \
        summary summary.dat \
        info info.dat \
        cpopvtime cpopvtime.agr normframe \
        repout {complex}_{time}_{epsilon} repfmt pdb \
        clusterout {complex} clusterfmt netcdf      
        ''')

    elif method == 'hieragglo':
        string = textwrap.dedent(f'''\
        trajin {traj} 1 last {interval}
        cluster {complex} \
        hieragglo epsilon {epsilon} complete epsilonplot {complex}_epsilon.dat \
        srmsd {mask},:{bs_mask} \
        out cnumvtime.dat \
        sil Sil \
        summary summary.dat \
        info info.dat \
        pairdist pairdist.dat \
        cpopvtime cpopvtime.agr normframe \
        repout rep repfmt pdb \
        singlerepout singlerep.nc singlerepfmt netcdf \
        avgout Avg avgfmt restart \
        clusterout {complex} clusterfmt netcdf
        ''')

    else:
        raise Exception(f'Method {method} not implemented')

    with open('cluster.ptraj', 'w') as f:
        f.write(string)

    run(f'cpptraj {parm} cluster.ptraj')

    os.chdir('../')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Cluster a trajectory by distance to ligand in binding site')

    # Required arguments
    parser.add_argument('-c','--complex', help='Name of complexes to run', required=True)
    parser.add_argument('-e', '--epsilon', type=float, help='Epsilon value to to clustering', required=True)

    # Optional arguments
    parser.add_argument('--exclude', default=None, nargs='+', type=int, help='Space separated list of resids to exclude', required=False)
    parser.add_argument('-f', '--folder_name', help='Name for folder', required=False)
    parser.add_argument('-i', '--interval', default=1, type=int, help='Interval to analyse trajectories', required=False)
    parser.add_argument('-m', '--mask', help='Mask to cluster around', required=False)
    parser.add_argument('--method', default='dbscan', help='Method used to do clustering, options are dbscan and hieragglo', required=False)
    parser.add_argument('-p', '--parm', help='Parameter file', required=False)
    parser.add_argument('-r','--repeat', default='', help='Repeat pattern, e.g. _2, _3', required=False)
    parser.add_argument('-t', '--time', default='equi', help='Time pattern to run second MD, etc, e.g. equi, equi2', required=False)
    parser.add_argument('--traj', help='Trajectory file', required=False)

    args = parser.parse_args()

    cluster_bs(complex=args.complex, epsilon=args.epsilon, exclude=args.exclude, folder_name=args.folder_name, interval=args.interval, mask=args.mask, method=args.method, parm=args.parm, repeat=args.repeat, time=args.time, traj=args.traj)

