import sys

sys.path.append('../')
from run import *

def test_run():
	print('line7')
	output = run('echo Hello')
	print('line9')
	print(output)

	assert output == 'Hello\n'

	try:
		print('line15')
		output = run('command_not_found')
	except:
		print('line18')
		print('Exception handled')
	print('line20')
	output = run('echo End')


def test_run2():

	run('bash test_run.sh')

test_run2()