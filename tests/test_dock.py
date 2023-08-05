import filecmp
import os
import sys

sys.path.append('../')
from dock import *

def test_get_first_mol2():

    get_first_mol2('input/6td3_1-1_dock.mol2', 'output/6td3_1-1_dock_1.mol2')

    assert filecmp.cmp('input/6td3_1-1_dock_1.mol2', 'output/6td3_1-1_dock_1.mol2')

def test_split_sdf():

    os.system('cp input/cycle_confs.sdf output/cycle_confs.sdf')

    n_files = split_sdf('output/cycle_confs.sdf')

    assert n_files == 11

    for i in range(n_files):

        assert filecmp.cmp(f'input/cycle_confs_{i}.sdf', f'output/cycle_confs_{i}.sdf')

def test_split_sdf_max1():

    os.system('cp input/cycle_confs.sdf output/cycle_confs.sdf')

    n_files = split_sdf('output/cycle_confs.sdf', max_n=1)

    print(n_files)

    assert n_files == 1

    for i in range(n_files):

        assert filecmp.cmp(f'input/cycle_confs_{i}.sdf', f'output/cycle_confs_{i}.sdf')

def test_split_sdf_max2():

    os.system('cp input/cycle_confs.sdf output/cycle_confs.sdf')

    n_files = split_sdf('output/cycle_confs.sdf', max_n=2)

    print(n_files)

    assert n_files == 2

    for i in range(n_files):

        assert filecmp.cmp(f'input/cycle_confs_{i}.sdf', f'output/cycle_confs_{i}.sdf')


def test_dock_conf():

    os.system('cp input/6td3_protein.pdb output/')

    dock_conf('1-1', 'input', 'output')


def test_parse_dock():

    score = parse_dock('input/1-1_confs_0_dock.out')

    assert score == -8.96

def test_write_gnina_output():

    score = parse_dock('input/1-1_confs_0_dock.out')

    write_gnina_output(score, 'output/1-1_confs_0_write.out')

    score2 = parse_dock('output/1-1_confs_0_write.out')

    assert score == score2

test_split_sdf_max2()

