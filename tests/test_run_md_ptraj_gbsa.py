import os

def test1():

	os.system('rm -rf output/*')

	os.system('cp input/l23.parm7 output/')
	os.system('cp input/l23.rst7 output/')

	os.chdir('output')

	os.system('python ../../run/run_md_ptraj_gbsa.py -c l23 --test YES')

	assert os.path.exists('gbsa/l23_gbsa.dat')

	os.chdir('../')


def test2():

	os.system('rm -rf output/*')

	os.system('cp input/l23.parm7 output/')
	os.system('cp input/l23.rst7 output/')

	os.chdir('output')

	os.system('python ../../run/run_md_ptraj_gbsa.py -c l23 --test YES -l L23')

	assert os.path.exists('gbsa/l23_gbsa.dat')

	os.chdir('../')

test2()