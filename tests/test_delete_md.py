import sys

sys.path.append('../')
from util.delete_md import *

def test_parse_count_frames_output1():
	# test for a simple string input

	n_frames = parse_count_frames_output('Read 500 frames and processed 500 frames.\n')

	assert n_frames == 500


def test_count_frames():
	# test count frames for a sample MD with 2 frames

	n_frames = count_frames('input/l23.parm7', 'input/l23_equi.nc')

	print(n_frames)

	assert n_frames == 2


def test_delete_mds():
	# test delete traj for sample MD with 2 frames

	os.system('cp input/l23.parm7 output/')
	os.system('cp input/l23_strip.parm7 output/')
	os.system('cp input/l23_equi.nc output/')
	os.system('cp input/l23_equi_cent.nc output/')
	os.system('cp input/l23_equi_cent_strip.nc output/')

	os.chdir('output')

	delete_mds('l23')

	assert os.path.isfile('l23_equi_cent.nc') is False
	assert os.path.isfile('l23_equi.nc') is False
	assert os.path.isfile('l23_equi_cent_strip.nc') is True

	os.chdir('../')


def test_main():
	# check main for sample MD with 2 frames

	os.system('cp input/l23.parm7 output/')
	os.system('cp input/l23_strip.parm7 output/')
	os.system('cp input/l23_equi.nc output/')
	os.system('cp input/l23_equi_cent.nc output/')
	os.system('cp input/l23_equi_cent_strip.nc output/')

	os.chdir('output')

	run(f'python ../../util/delete_md.py -i l23')	

	assert os.path.isfile('l23_equi_cent.nc') is False
	assert os.path.isfile('l23_equi.nc') is False
	assert os.path.isfile('l23_equi_cent_strip.nc') is True

	os.chdir('../')


def test_main_repeat():
	# check main for sample MD with 2 frames and repeat pattern

	repeat = '_2'

	os.system('cp input/l23.parm7 output/')
	os.system('cp input/l23_strip.parm7 output/')
	os.system(f'cp input/l23_equi.nc output/l23{repeat}_equi.nc')
	os.system(f'cp input/l23_equi_cent.nc output/l23{repeat}_equi_cent.nc')
	os.system(f'cp input/l23_equi_cent_strip.nc output/l23{repeat}_equi_cent_strip.nc')

	os.chdir('output')

	run(f'python ../../util/delete_md.py -i l23 -r {repeat}')	

	assert os.path.isfile(f'l23{repeat}_equi_cent.nc') is False
	assert os.path.isfile(f'l23{repeat}_equi.nc') is False
	assert os.path.isfile(f'l23{repeat}_equi_cent_strip.nc') is True

	os.chdir('../')