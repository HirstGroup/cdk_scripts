import os

def test1():
	# test for only part 2

	os.system('rm -rf output/*')

	os.system('cp input/l23.parm7 output/')
	os.system('cp input/l23_equi.nc output/l23_equi2.nc')

	os.chdir('output')

	os.system('python ../../run/run_ptraj_gbsa.py -c l23 -fp 2 -lp 2 --test YES')

	os.chdir('../')

	assert os.path.exists('output/gbsa2/l23_gbsa2.dat') is True


def test2():
	# test for parts 2 to 3 and repeat _2

	os.system('rm -rf output/*')

	os.system('cp input/l23.parm7 output/')
	os.system('cp input/l23_equi.nc output/l23_2_equi2.nc')
	os.system('cp input/l23_equi.nc output/l23_2_equi3.nc')

	os.chdir('output')

	os.system('python ../../run/run_ptraj_gbsa.py -c l23 -fp 2 -lp 3 -r _2 --test YES')

	os.chdir('../')

	assert os.path.exists('output/gbsa2_2/l23_gbsa2_2.dat') is True
	assert os.path.exists('output/gbsa3_2/l23_gbsa3_2.dat') is True

test2()