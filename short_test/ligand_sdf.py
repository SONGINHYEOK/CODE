import openbabel
from openbabel import openbabel
obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("pdb", "mol2")

mol = openbabel.OBMol()
obConversion.ReadFile(mol, "/Users/song-inhyeok/Documents/PDB/ligand/lig_6v43_H_LUZ_14.pdb")   # Open Babel will uncompress automatically


print(mol.NumAtoms())
print(mol.NumBonds())
print(mol.NumResidues())

obConversion.WriteFile(mol, '/Users/song-inhyeok/Documents/PDB/ligand/lig_6yp7_H_LUT_14.mol')




