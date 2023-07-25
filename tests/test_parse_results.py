import os
import sys

sys.path.append('../')
from parse_results import *

def test_parse_gbsa():

	check = {'vdwaals': -38.7039, 'eel': -186.3014, 'egb': 206.2819, 'esurf': -4.554, 'delta_g_gas': -225.0053, 'delta_g_solv': 201.7279, 'delta_total': -23.2774}

	result = parse_gbsa('input/a01_gbsa.dat')

	print(result)

	assert result == check

test_parse_gbsa()

def test_parse_gbsa_df():

	os.chdir('input')

	df = pd.DataFrame(['a01', 'a06'], columns=['ligname'])

	print(df)

	df = parse_gbsa_df(df)

	print(df)

	assert list(df['gbsa_delta_total']) == [-23.2774, -23.2774]

	os.chdir('../')

test_parse_gbsa_df()