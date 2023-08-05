import filecmp
import os
import sys

sys.path.append('../')
from ptraj import *

from compare_files import compare_files


def test_center_strip():

	os.system('cp input/a01.parm7 output/')

	os.system('cp input/a01_equi.nc output/')

	os.chdir('output')

	os.system('rm -rf a01_strip.parm7 a01_equi_cent_strip.nc')

	center_strip('a01')

	os.chdir('../')

	compare_files('input/a01_strip.parm7', 'output/a01_strip.parm7', 1)

	assert filecmp.cmp('input/a01_equi_cent_strip.nc', 'output/a01_equi_cent_strip.nc')


def test_ptraj_main_delete():
	# test ptraj main function with delete option for sample MD with 2 frames

	os.system('cp input/l23.parm7 output/')
	os.system('cp input/l23_equi.nc output/')
	os.system('rm output/l23_equi_cent.nc')
	os.system('rm output/l23_equi_cent_strip.nc')

	os.chdir('output')

	os.system('python ../../ptraj.py -i l23 --delete')	

	assert os.path.isfile('l23_equi_cent.nc') is False
	assert os.path.isfile('l23_equi.nc') is False
	assert os.path.isfile('l23_equi_cent_strip.nc') is True

	os.chdir('../')

test_ptraj_main_delete()