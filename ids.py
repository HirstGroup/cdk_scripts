infile = open('CDK12_CDK9_Original_Data14.smi')

outfile = open('CDK12_CDK9_Original_Data14_ids.smi', 'w')

n = 0

for line in infile:
    n += 1
    smi = line.split()[0].strip()
    
    outfile.write(f'{smi} {n}\n')
