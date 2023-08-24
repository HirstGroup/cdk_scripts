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


def get_bs_mask(frames, mask, parm, traj, interval=1):
    """
    Get mask of residues in binding site
    """

    if not os.path.exists('mask'):
        os.makedirs('mask')

    os.chdir('mask')

    string = textwrap.dedent(f'''\
    trajin ../{traj} 1 last {interval}
    mask :{mask}<:4.0 maskpdb resid_4a.pdb
    ''')

    with open('mask.ptraj', 'w') as f:
        f.write(string)

    run(f'cpptraj ../{parm} mask.ptraj', print_cmd=True)

    resid_set = set()

    for frame in range(1,frames+1):
        df = np.genfromtxt(f'resid_4a.pdb.{frames}', names=None, skip_header=1, usecols=4, skip_footer=1, dtype=int)
        resid_set.update(set(df))
    
    amber_range = get_amber_range(resid_set, 'resid_vmd.txt')

    os.chdir('../')

    return amber_range


def cluster_bs(complex, epsilon=1.0, interval=1, mask=None, method='dbscan', repeat='', time='equi'):
    """
    Cluster trajectory by distance to ligand in binding site
    """

    if mask is None:
        mask = complex.upper()

    parm = f'../{complex}_strip.parm7'
    traj = f'../{complex}{repeat}_{time}_cent_strip.nc'

    folder = f'cluster{repeat}_{time}_{interval}_{epsilon}'

    if not os.path.exists(folder):
        os.makedirs(folder)

    os.chdir(folder)

    frames = count_frames(parm, traj, verbose=True)

    mask = get_bs_mask(frames, mask, parm, traj, interval=interval)

    if method == 'dbscan':
        string = textwrap.dedent(f'''\
        trajin {traj} 1 last {interval}
        cluster {complex} \
        dbscan minpoints 2 epsilon {epsilon} sievetoframe \
        rms :{mask} nofit \
        sieve 1 random \
        pairdist pairdist.dat \
        out cnumvtime.dat \
        sil Sil \
        summary summary.dat \
        info info.dat \
        cpopvtime cpopvtime.agr normframe \
        repout {complex}_{time}_{epsilon} repfmt pdb       
        ''')

    elif method == 'hieragglo':
        string = textwrap.dedent(f'''\
        trajin {traj} 1 last {interval}
        cluster {complex} \
        hieragglo epsilon {epsilon} complete epsilonplot {complex}_epsilon.dat \
        srmsd * \
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
    parser.add_argument('-i', '--interval', default=1, type=int, help='Interval to analyse trajectories', required=False)
    parser.add_argument('-m', '--mask', help='Mask to cluster around', required=False)
    parser.add_argument('--method', default='dbscan', help='Method used to do clustering', required=False)
    parser.add_argument('-r','--repeat', default='', help='Repeat pattern, e.g. _2, _3', required=False)
    parser.add_argument('-t', '--time', default='equi', help='Time pattern to run second MD, etc, e.g. equi, equi2', required=False)

    args = parser.parse_args()

    cluster_bs(complex=args.complex, epsilon=args.epsilon, interval=args.interval, mask=args.mask, method=args.method, repeat=args.repeat, time=args.time)

