import filecmp
import sys

sys.path.append('../')
from ptraj_scripts.cluster_bs import *


def test_get_range_tuple():

	resid_list = [1,2,3,6,7,10]

	range_tuple = list(get_range_tuple(resid_list))

	assert range_tuple == [(1, 3), (6, 7), (10, 10)]	


def test_get_amber_ranges():

	resid_list = [1,2,3,6,7,10]

	amber_range = get_amber_range(resid_list, 'output/range_file.txt')

	print(amber_range)

	assert amber_range == '1-3,6-7,10,'

	assert filecmp.cmp('input/range_file.txt', 'output/range_file.txt') is True


def test_get_bs_mask():

	os.system('rm -rf output/*')
	os.system('cp input/l23_strip.parm7 output/')
	os.system('cp input/l23_equi_cent_strip.nc output/')

	os.chdir('output/')

	mask = get_bs_mask(2, ':L23', 'l23_strip.parm7', 'l23_equi_cent_strip.nc', interval=1)

	print(mask)

	os.chdir('../')

	assert mask == '57-58,61,65,76,79-80,83-84,93,96-97,100-101,104,120,124,151,'


def test_cluster_bs():

	os.system('rm -rf output/*')
	os.system('cp input/l23_strip.parm7 output/')
	os.system('cp input/l23_equi_cent_strip_10.nc output/l23_equi_cent_strip.nc')

	os.chdir('output')

	cluster_bs('l23')

	os.chdir('../')

test_get_range_tuple()

test_get_amber_ranges()

test_get_bs_mask()

test_cluster_bs()
