from openbabel import openbabel
from openbabel import pybel


obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("sdf","pdb")

mol = openbabel.OBMol()

name  = 'MD_tg_mol'
     
root = "/Users/song-inhyeok/Documents/PDB/ligand/"

obConversion.ReadFile(mol, root+"1V4S_A_GLC_500_1.sdf") # Open Babel will uncompress automatically

obConversion.WriteFile(mol,name+".pdb")

