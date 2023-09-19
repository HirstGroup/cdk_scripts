from openeye import oechem, oeff

ifs = oechem.oemolistream()
ifs.open("input/a01.mol2")

for mol in ifs.GetOEMols():

	# Set up the MMFF94 force field
	forcefield = oeff.OEForceFieldType_MMFF94()

	# Calculate the energy
	energy = forcefield.CalcEnergy(mol)
	print("Molecular Energy:", energy)