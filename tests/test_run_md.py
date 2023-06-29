import filecmp
import sys

sys.path.append('../')
from run_md import *

from pandas.testing import assert_frame_equal

def test_number_to_resname():

	assert number_to_resname(1) == 'A01'
	assert number_to_resname(2) == 'A02'
	assert number_to_resname(10) == 'A10'
	assert number_to_resname(100) == 'B00'
	assert number_to_resname(200) == 'C00'
	assert number_to_resname(201) == 'C01'
	assert number_to_resname(299) == 'C99'


def test_make_resname():

	d = {'Row':[1,2,3], 'stereoisomers_list':['1&2','1&2&3','1']}

	df = pd.DataFrame(data=d)

	print(df)

	df = make_resname(df)

	print(df)

	d2 = {'Row':[1,2,3], 'stereoisomers_list':['1&2','1&2&3','1'], 'resname_list':['A01&A02', 'A03&A04&A05', 'A06']}

	df2 = pd.DataFrame(data=d2)

	assert_frame_equal(df, df2)



test_make_resname()
