#!/usr/bin/env python

import argparse
import os
import pandas as pd


def get_results(row, aux):
    """
    Get average results from each ligname into a dictionary

    Parameters
    ----------
    row: pandas df
    aux: pandas df with data to extract

    Returns
    -------
    results: df_out
        pandas dataframe with extracted data
    """

    results = dict()

    resname_list = row['resname_list'].split('&')

    resname_list_lower = [i.lower() for i in resname_list]

    aux = aux[aux['ligname'].isin(resname_list_lower)]

    aux = aux.drop('ligname', axis=1)

    aux = aux.mean().to_dict()

    return aux


def main(input, aux, output):
    """
    Combine results for each ligname into an experimental result

    Parameters
    ----------
    input: str
        input file name csv file separated by semicolon
    aux: str
        auxiliary file name csv file separated by semicolon where to get results from
    output: str
        output file file name csv file separated by semicolon where results will be saved

    Returns
    -------
    None (creates output file)
    """

    if input == output:
        os.system(f'cp {input} {input}.bk')

    df = pd.read_csv(input, sep=';')

    aux = pd.read_csv(aux, sep=';')

    df['results'] = ''

    for index, row in df.iterrows():
        df.at[index,'results'] = get_results(row, aux)

    df = df.join(pd.json_normalize(df['results'])).drop('results', axis=1)

    df.to_csv(output, sep=';', index=False)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Combine results for each ligname into an experimental result')

    # Required arguments
    parser.add_argument('-i','--input', help='Input file name, csv file separated by semicolon',required=True)
    parser.add_argument('-a','--aux', help='Auxiliary file name, csv file separated by semicolon with results to add to input data',required=True)
    parser.add_argument('-o','--output', help='Output file name',required=True)

    args = parser.parse_args()

    if args.input == args.output:
        os.system(f'cp {args.input} {args.input}.bk')    

    main(args.input, args.aux, args.output)





