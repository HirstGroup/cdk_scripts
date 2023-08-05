import filecmp
import os
import sys

from pandas.testing import assert_frame_equal

sys.path.append('../')
from add_results import *

def test_get_results():


	df = pd.DataFrame(['A16&A17'], columns=['resname_list'])

	print(df)

	aux = pd.read_csv('input/ligname-data.csv', sep=';')

	print(aux)

	df['results'] = ''

	for index, row in df.iterrows():
		df.at[index,'results'] = get_results(row, aux)

	print(df.to_string())

	check = pd.DataFrame(zip(['A16&A17'], [{'vdwaals': -34.515, 'eel': -35.781, 'egb': 37.63915, 'esurf': -4.30715, 'delta_g_gas': -70.29599999999999, 'delta_g_solv': 33.332, 'delta_total': -36.964}]), columns=['resname_list', 'results'])

	print(check)

	assert_frame_equal(df, check)


def test_main():

	main('input/Nottingham_data15.csv', 'input/ligname-data.csv', 'output/Nottingham_data16.csv')

	assert filecmp.cmp('input/Nottingham_data16.csv', 'output/Nottingham_data16.csv')

test_main()

