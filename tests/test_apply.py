import pandas as pd


def mass_times_date(row, multiple=1):

	return row['mass'] * row['date'] * multiple


a = 'mol'
b = 'mass'
c = 'date'


a = [{a:1, b:2, c:3}, {a:4, b:5, c:6}]

df = pd.DataFrame(a)

print(df)

df['func'] = df.apply(mass_times_date, axis=1)

df['func2'] = df.apply(mass_times_date, axis=1, multiple=2)

print(df)

