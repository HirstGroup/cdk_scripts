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


def test_ignore_errors():
	# execute run ignoring errors, i.e. don't crash

	# this command will crash
	try:
		run('command_not_found')
	except:
		print('Command crashed as expected')

	# this command will crash but will be ignored
	run('command_not_found', ignore_errors=True)
	print('command crashed but ignored')
