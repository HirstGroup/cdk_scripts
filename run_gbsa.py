from gbsa import *

def create_resp1_file_row(row):

	for x, smi in enumerate(row['stereoisomers_list'].split('&')):

		ligname = '%s-%s' %(row['Row'], x+1)

		name = row['Row']

		infile = ligname + '.mol2'

		outfile = 'resp/' + ligname + '_opt.gau'

		charge = get_charge(infile)

		create_resp1_file(infile, outfile, charge, cpu=10)


df = pd.read_csv('Nottingham_data12.csv', sep=';')

df = df.loc[df['Covalent'] == False]

df.dropna(inplace=True, subset=['CDK12 Mean IC50 (uM)'])

df.sort_values(by='Row', inplace=True)

df = df.head(n=20)

df.apply(create_resp1_file_row, axis=1)

