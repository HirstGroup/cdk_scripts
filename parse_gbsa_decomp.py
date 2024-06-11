import pandas as pd


def parse_gbsa_decomp(input, output=None):
	"""
	Parse GBSA decomp ASCII output to csv file and pandas df

	Parameters
	----------
	input : str
		Input file name
	output : str
		Output CSV file name

	Returns
	-------
	df : Pandas df
		Pandas df with GBSA decomp data
	"""

	with open(input) as f:

		lines = f.readlines()

	for x, line in enumerate(lines):
		if 'Total Energy Decomposition' in line:
			first = x + 3

		if 'Sidechain Energy Decomposition' in line:
			last = x - 1

	lines = lines[first:last]

	resname_list = []
	resid_list = []
	total_list = []

	for line in lines:

		resname = line.split()[0]
		resid = int(line.split()[1])
		total = float(line.split()[27])

		resname_list.append(resname)
		resid_list.append(resid)
		total_list.append(total)

	df = pd.DataFrame(list(zip(resname_list, resid_list, total_list)), columns =['Resname', 'Resid', 'Total'])

	if output is not None:
		df.to_csv(output, sep=',', index=False)

	return df


def average_decomp(input_list, output=None):
	"""
	Take averages of GBSA decomp outputs

	Parameters
	----------
	input : list of str
		List of input file names
	output : str
		Output CSV file name	

	Returns
	-------
	df : Pandas df
		Pandas df with average GBSA decomp data	
	"""

	df = parse_gbsa_decomp(input=input_list[0])
	df = df.rename(columns={'Total': 'Total_0'})

	column_names = ['Total_0']

	for x, f in enumerate(input_list[1:]):

		df2 = parse_gbsa_decomp(f)

		assert list(df['Resname']) == list(df2['Resname'])
		assert list(df['Resid']) == list(df2['Resid'])

		name = f'Total_{x+1}'
		column_names.append(name)

		df[name] = df2['Total']

	df['Total'] = df[column_names].mean(axis=1)
	df['Total_std'] = df[column_names].std(axis=1)

	if output is not None:
		df.to_csv(output, sep=',', index=False)

	return df



