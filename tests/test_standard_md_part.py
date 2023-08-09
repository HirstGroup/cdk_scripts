import sys

sys.path.append('../')
from standard_md_part import *

def test1():
	# test with part = 2 for l23

	os.system('rm output/*')
	os.system('cp input/l23.parm7 output/')
	os.system('cp input/l23_heat2.rst7 output/')

	standard_md_part(complex='l23', part=2, cd='output', test='YES')

	with open('output/l23_equi2.out') as f:
		assert 'Total wall time' in f.read()

def test2():
	# test with part = 3 and repeat = _2 for l23

	os.system('rm output/*')
	os.system('cp input/l23.parm7 output/')
	os.system('cp input/l23_equi.rst7 output/l23_2_equi2.rst7')

	standard_md_part(complex='l23', part=3, cd='output', repeat='_2', test='YES')

	with open('output/l23_2_equi3.out') as f:
		assert 'Total wall time' in f.read()

def test3():
	# test with part = 3 and repeat = NO for l23

	os.system('rm output/*')
	os.system('cp input/l23.parm7 output/')
	os.system('cp input/l23_equi.rst7 output/l23_equi2.rst7')

	standard_md_part(complex='l23', part=3, cd='output', repeat='NO', test='YES')

	with open('output/l23_equi3.out') as f:
		assert 'Total wall time' in f.read()

test3()