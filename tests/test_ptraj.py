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

test_center_strip()