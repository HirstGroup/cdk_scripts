import filecmp
import sys

sys.path.append('../')
from util.check_run import *

def test1():

	os.system('rm -rf output/*')
	os.system('cp input/l23.parm7 output/')
	os.system('cp input/l23_strip.parm7 output/')
	os.system('cp input/l23_equi.nc output/')
	os.system('cp input/l23_equi_cent.nc output/')
	os.system('cp input/l23_equi_cent_strip.nc output/')
	os.system('mkdir output/gbsa/')
	os.system('cp input/l23_gbsa2.dat output/gbsa/l23_gbsa.dat')

	os.chdir('output')

	os.system('python ../../util/check_run.py -i l23 -f 2 -o test_check_run1.txt')

	os.chdir('../')

	assert filecmp.cmp('input/test_check_run1.txt', 'output/test_check_run1.txt') is True


def test2():

	os.system('rm -rf output/*')
	os.system('cp input/l23.parm7 output/')
	os.system('cp input/l23_strip.parm7 output/')
	
	os.system('cp input/l23_equi.nc output/l23_2_equi.nc')
	os.system('cp input/l23_equi_cent.nc output/l23_2_equi_cent.nc')
	os.system('cp input/l23_equi_cent_strip.nc output/l23_2_equi_cent_strip.nc')
	os.system('mkdir output/gbsa_2/')
	os.system('cp input/l23_gbsa2.dat output/gbsa_2/l23_gbsa_2.dat')

	os.system('cp input/l23_equi.nc output/l23_2_equi2.nc')
	os.system('cp input/l23_equi_cent.nc output/l23_2_equi2_cent.nc')
	os.system('cp input/l23_equi_cent_strip.nc output/l23_2_equi2_cent_strip.nc')
	os.system('mkdir output/gbsa2_2/')
	os.system('cp input/l23_gbsa2.dat output/gbsa2_2/l23_gbsa2_2.dat')

	os.chdir('output')

	os.system('python ../../util/check_run.py -i l23 -f 2 -p 1 2 -o test_check_run2.txt -r _2 --verbose')

	os.chdir('../')

	assert filecmp.cmp('input/test_check_run2.txt', 'output/test_check_run2.txt') is True