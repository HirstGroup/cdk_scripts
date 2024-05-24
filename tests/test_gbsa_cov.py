import filecmp
import os
import sys

sys.path.append('../')
from gbsa_scripts.gbsa_cov import *

def test_determine_strip_atoms():

	os.chdir('output')

	os.system('cp ../input/ads081.parm7 ../input/ads081.rst7 .')

	strip_atoms = determine_strip_atoms(complex='ads081', cov_indeces=[4,62,3,2,61,0,1,64,63,65], ligandmask='F01')

	assert strip_atoms == [5152, 5153, 5154, 5155, 5156, 5157, 5158, 5159, 5160, 5161, 5162, 5163, 5164, 5165, 5166, 5167, 5168, 5169, 5170, 5171, 5172, 5173, 5174, 5175, 5176, 5177, 5178, 5179, 5180, 5181, 5182, 5183, 5184, 5185, 5186, 5187, 5188, 5189, 5190, 5191, 5192, 5193, 5194, 5195, 5196, 5197, 5198, 5199, 5200, 5201, 5202, 5203, 5204, 5205, 5206, 5207]

