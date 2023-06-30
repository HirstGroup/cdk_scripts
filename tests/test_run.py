import sys

sys.path.append('../')
from run import *

def test_run():

	output = run('echo Hello')

	assert output == 'Hello\n'

	try:
		output = run('command_not_found')
	except:
		print('Exception handled')

	output = run('echo End')

test_run()