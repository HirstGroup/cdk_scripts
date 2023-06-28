#/usr/bin/env python

import argparse
import operator

from gbsa import *


def create_resp1_file_row(row):

    for x, ligname in enumerate(row['resname_list'].split('&')):

        ligname = ligname.lower()

        infile = 'dock_conf/' + ligname + '_dock_best.sdf'

        outfile = 'gbsa/resp/' + ligname + '_opt.gau'

        charge = get_charge(infile)

        create_resp1_file(infile, outfile, charge, cpu=10)


def create_resp1_file_row_simple(row):

    ligname = row['ligname']

    infile = 'dock_conf/' + ligname + '_dock_best.sdf'

    outfile = 'resp/' + ligname + '_opt.gau'

    charge = get_charge(infile)

    create_resp1_file(infile, outfile, charge, cpu=10)


def check_resp1_output_row(row):

    for x, smi in enumerate(row['stereoisomers_list'].split('&')):

        ligname = '%s-%s' %(row['Row'], x+1)

        infile = 'gbsa/resp/' + ligname + '_opt.log'

        inchikey = get_inchikey(ligname + '.mol2')

        check = check_resp1_output(infile, inchikey)

        return check


def get_charge_row(row):

    charge_list = []

    for x, smi in enumerate(row['stereoisomers_list'].split('&')):

        ligname = '%s-%s' %(row['Row'], x+1)

        infile = ligname + '.mol2'

        outfile = 'resp/' + ligname + '_opt.gau'

        charge_list.append( get_charge(infile) )

    assert len(set(charge_list)) == 1

    return charge_list[0]


def create_resp2_file_row(row):

    for x, smi in enumerate(row['stereoisomers_list'].split('&')):

        ligname = '%s-%s' %(row['Row'], x+1)

        infile = 'gbsa/resp/' + ligname + '_opt.log'

        outfile = 'gbsa/resp/' + ligname + '_esp.gau'

        create_resp2_file(infile, outfile, row['charge'], cpu=8)


def create_resp3_file_row(row):

    for x, smi in enumerate(row['stereoisomers_list'].split('&')):

        ligname = '%s-%s' %(row['Row'], x+1)

        infile = 'gbsa/resp/' + ligname + '_esp.log'

        outfile = 'gbsa/resp/' + ligname + '_resp.mol2'

        auxfile = ligname + '.mol2'

        inchikey = get_inchikey(ligname + '.mol2')

        resname = inchikey[0:3]

        create_resp3_file(infile, outfile, auxfile, resname)


def make_resname(df):

    df.sort_values(by='Row', inplace=True)

    print(df)

    count = 0

    resname_list = []

    for index, row in df.iterrows():

        resname_row = []

        for x, smi in enumerate(row['stereoisomers_list'].split('&')):

            count +=1

            print(count)

            resname_row.append(number_to_resname(count))

        print('&'.join(resname_row))

        resname_list.append('&'.join(resname_row))


    df['resname_list'] = resname_list

    return df


def number_to_resname(number):

    resname = str(number).zfill(3)

    d = {'0':'A', '1':'B', '2':'C', '3':'D', '4':'E', '5':'F', '6':'G'}

    resname = d[resname[0]] + resname[1:3]

    return resname


def rename_ligands(row):

    for x, smi in enumerate(row['stereoisomers_list'].split('&')):

        name = '%s-%s' %(row['Row'], x+1)

        newname = row['resname_list'].split('&')[x].lower()

        cmd = f'rename "s/{name}/{newname.lower()}/" *_{name}_*.*'

        #cmd = f'mv {name}-3d.sdf {newname}-3d.sdf'

        print(cmd)

        os.system(cmd)


def print_ligand(row):

    print(row['ligname'])


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run gbsa.py')

    # Required arguments
    parser.add_argument('-f','--function', help='Function name to run',required=True)
    parser.add_argument('-i','--input', help='Input file name, csv file separated by semicolon if using pandas option, otherwise space separated list of ligand names',required=True)

    # Optional arguments
    parser.add_argument('-o','--output', help='Output file name',required=False)
    parser.add_argument('--pandas', action='store_true', help='Apply Pandas function row by row', required=False)

    args = parser.parse_args()

    if args.input == args.output:
        os.system(f'cp {args.input} {args.input}.bk')

    if args.pandas:
        df = pd.read_csv(args.input, sep=';')

        #df = df.loc[df['Covalent'] == False]

        #df.dropna(inplace=True, subset=['CDK12 Mean IC50 (uM)'])

        #df.sort_values(by='Row', inplace=True)

        #df = df.head(n=20)

        df[args.function] = df.apply(eval(args.function), axis=1)

        # remove added column if function has None as output
        col = df[args.function].tolist()

        if len(set(col)) == 1 and col[0] is None:
            df.drop([args.function], axis=1, inplace=True)

        if args.output is not None:
            df.to_csv(args.output, sep=';', index=False)

    else:

        df = pd.DataFrame(args.input.split(), columns=['ligname'])

        print(df)

        df.apply(eval(args.function), axis=1)