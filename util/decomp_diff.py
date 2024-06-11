import pandas as pd
import sys

complex = sys.argv[1]

df = pd.read_csv(f'./redock_left/{complex}/{complex}_gbsa_decomp_avg.csv')

df = df[2:330]
df = df.reset_index(drop=True)

df2 = pd.read_csv(f'./redock_right/{complex}/{complex}_gbsa_decomp_avg.csv')

df2 = df2[:328]
df2 = df2.reset_index(drop=True)

assert list(df['Resname']) == list(df2['Resname'])

df['Resname_sum'] = df['Resname'] + df2['Resname']

df['Resname_left'] = df['Resname']
df['Resname_right'] = df2['Resname']

df['Resid_left'] = df['Resid']
df['Resid_right'] = df2['Resid']

assert list(df['Resname_left']) == list(df['Resname_right'])

df['Total_left'] = df['Total']
df['Total_right'] = df2['Total']

df['Total_left_std'] = df['Total_std']
df['Total_right_std'] = df2['Total_std']

df['diff'] = df['Total_right'] - df['Total_left']
df['diff_abs'] = df['diff'].abs()

df = df.sort_values(by=['diff_abs'], ascending=False)

df.to_csv(f'{complex}_diff.csv', sep=',', index=False)

print(df)