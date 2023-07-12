import argparse 

parser = argparse.ArgumentParser(description='Grep inchi')
parser.add_argument('-i','--input', help='Input File',required=True)
parser.add_argument('-a','--auxiliary_inchi_file', help='Auxiliary Inchi file containing Inchi patterns to grep',required=True)

args = parser.parse_args()

infile = open(args.input)

auxfile = open(args.auxiliary_inchi_file)

inchi_list = []

for line in auxfile:
    inchi = line.split()[0]
    
    inchi_list.append(inchi)


for line in infile:
    for inchi in inchi_list:
        if inchi in line:
            print(line.strip())


