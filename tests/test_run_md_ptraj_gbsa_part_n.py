import os

def test1():
	# Test for part 2 and repeat _3

	os.system('rm -rf output/*')
	os.system('cp input/l23.parm7 output/')
	os.system('cp input/l23_heat2.rst7 output/l23_3_equi.rst7')

	os.chdir('output')

	os.system('python ../../run/run_md_ptraj_gbsa_part_n.py -c l23 -fp 2 -lp 2 -r _3 --sequential --test YES')

	i = 2
	for file_name in f'l23_3_equi{i}.out l23_3_equi{i}.rst7 l23_3_equi{i}_cent_strip.nc gbsa{i}_3/l23_gbsa{i}_3.dat'.split():
		assert os.path.exists(file_name)

	os.chdir('../')


def test2():
	# Test for parts 2 to 10 and repeat _3

	os.system('rm -rf output/*')
	os.system('cp input/l23.parm7 output/')
	os.system('cp input/l23_heat2.rst7 output/l23_3_equi.rst7')

	os.chdir('output')

	os.system('python ../../run/run_md_ptraj_gbsa_part_n.py -c l23 -fp 2 -lp 10 -r _3 --sequential --test YES')

	for i in range(2,11):
		for file_name in f'l23_3_equi{i}.out l23_3_equi{i}.rst7 l23_3_equi{i}_cent_strip.nc gbsa{i}_3/l23_gbsa{i}_3.dat'.split():
			assert os.path.exists(file_name)

	os.chdir('../')


test2()
