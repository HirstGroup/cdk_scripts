import pandas as pd

df = pd.read_csv('Nottingham_data10.csv', sep=';')

df.sort_values(by='smi_flat')

d = dict()

for index, row in df.iterrows():

	if row['smi_flat'] in d:
		d[row['smi_flat']].append(row['Row'])

	else:
		d[row['smi_flat']] = [row['Row']]


outfile = open('figure.tex', 'w')

outfile.write('''
\\documentclass{article}
\\usepackage[utf8]{inputenc}
\\usepackage{graphicx}
\\usepackage{tikz}
\\usepackage[a4paper, margin=1.5cm]{geometry}
\\usepackage{pdfpages}

\\begin{document}

''')

a = []

for key, val in d.items():

	if len(val) == 1:
		continue

	outfile.write('\\newline\n')

	for i in val:
		outfile.write('%s - ' %i)

		a.append(i)

	outfile.write('\n')

	outfile.write('\\newline\n')

	for i in val:

		outfile.write('\includegraphics[scale=0.5]{fig/%s_crop.pdf}\n' %i)

outfile.write('\end{document}\n')


print(a)

df2 = df[~df['Row'].isin(a)]

df2.to_csv('Nottingham_data10_notrepeat.csv', sep=';', index=False)