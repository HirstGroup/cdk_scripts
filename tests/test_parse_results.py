import os
import sys

sys.path.append('../')
from parse_results import *

def test_parse_gbsa():

	results = {'vdwaals': -38.7039, 'eel': -186.3014, 'egb': 206.2819, 'esurf': -4.554, 'delta_g_gas': -225.0053, 'delta_g_solv': 201.7279, 'delta_total': -23.2774}

	for prop in ['vdwaals', 'eel', 'egb', 'esurf', 'delta_g_gas', 'delta_g_solv', 'delta_total']:
		result = parse_gbsa('input/a01_gbsa.dat', prop)
		assert results[prop] == result


def test_parse_gbsa_df():

	os.chdir('input')

	df = pd.DataFrame(['a01', 'a06'], columns=['ligname'])

	print(df)

	df = parse_gbsa_df(df)

	print(df)

	assert list(df['delta_total']) == [-23.2774, -23.2774]

	os.chdir('../')

test_parse_gbsa_df()