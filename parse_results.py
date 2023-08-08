#!/usr/bin/env python

import argparse
import os
import pandas as pd


def parse_gbsa(infile):
    """
    Parse GBSA output data

    Parameters
    ----------
    infile: str
        name of input file to parse
    prop: str
        property to return

    Returns
    -------
    results: float
        parsed result of prop
    """

    results = dict()

    sel = False

    with open(infile) as f:
        for line in f:
            if sel:
                if 'DELTA TOTAL' in line:
                    results['delta_total'] = float(line.split()[2])
                if 'VDWAALS' in line:
                    results['vdwaals'] = float(line.split()[1])
                if 'EEL' in line:
                    results['eel'] = float(line.split()[1])
                if 'EGB' in line:
                    results['egb'] = float(line.split()[1])
                if 'ESURF' in line:
                    results['esurf'] = float(line.split()[1])
                if 'DELTA G gas' in line:
                    results['delta_g_gas'] = float(line.split()[3])
                if 'DELTA G solv' in line:
                    results['delta_g_solv'] = float(line.split()[3])

            if 'Differences' in line:
                sel = True

    return results


def parse_gbsa_row(row, repeat='', interval=''):

    ligname = row['ligname']

    try:
        results = parse_gbsa(f'md/{ligname}/gbsa{repeat}/{ligname}_gbsa{repeat}{interval}.dat')

        new_results = {}

        for key, val in results.items():
            new_results[f'gbsa{repeat}{interval}_{key}'] = val

        return new_results
    except:
        print(f'No results for {ligname}')
        return None

def parse_gbsa_df(df, repeat='', interval=''):
    """
    Parse GBSA output data row by row.
    If new data columns exist in original df, original df columns will be removed.

    Parameters:
    df: pandas df

    Returns
    -------
    df: pandas df
    """

    df2 = pd.DataFrame()

    df2[f'gbsa_results{repeat}{interval}'] = df.apply(parse_gbsa_row, repeat=repeat, interval=interval, axis=1)

    df2 = df2.join(pd.json_normalize(df2[f'gbsa_results{repeat}{interval}'])).drop(f'gbsa_results{repeat}{interval}', axis=1)

    for col in df2.columns:
        if col in df.columns:
            df.drop(col, axis=1, inplace=True)

    df = df.join(df2)

    return df


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Parse results')

    # Required arguments
    parser.add_argument('-i','--input', help='Input file name, csv file separated by semicolon',required=True)

    # Optional arguments
    parser.add_argument('-o','--output', help='Output file name',required=False)
    parser.add_argument('-r','--repeat', default='', help='Repeat pattern to parse gbsa, e.g. _2, _3 etc.',required=False)
    parser.add_argument('--interval', default='', help='Interval pattern to parse gbsa, e.g. _2, _3 etc.',required=False)

    args = parser.parse_args()

    if args.input == args.output:
        os.system(f'cp {args.input} {args.input}.bk')

    df = pd.read_csv(args.input, sep=';')

    #df = df.loc[df['Covalent'] == False]

    #df.dropna(inplace=True, subset=['CDK12 Mean IC50 (uM)'])

    #df.head(n=15).apply(generate_mol2_row, axis=1)

    df = parse_gbsa_df(df, args.repeat, args.interval)

    print(df)

    if args.output is not None:
        df.to_csv(args.output, sep=';', index=False)
