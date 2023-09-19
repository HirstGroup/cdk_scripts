import filecmp
import os
import sys

sys.path.append('../')
from dock_scripts.dock_openeye import *

def test1():
	# Test full script 

	os.system('cp input/1-1.sdf output/')
	os.system('cp input/6td3_protein.oeb output/')
	
	os.chdir('output')

	oe_dock('1-1.sdf', '6td3_protein.oeb')

	os.chdir('../')

test1()