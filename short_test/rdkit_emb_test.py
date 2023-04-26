from rdkit import Chem, RDConfig
from rdkit.Chem import Draw, AllChem, rdMolAlign, MolStandardize



mol = Chem.MolFromPDBFile("/Users/song-inhyeok/Documents/PDB/ligand/lig_6v43_H_LUZ_14.pdb")
mol2 = Chem.AddHs(mol)
AllChem.ConstrainedEmbed(mol2, mol)
mol_3d = Chem.MolToMolBlock(mol)





def ligand_emb(mol):
    mol2 = Chem.AddHs(mol)
    AllChem.ConstrainedEmbed(mol2, mol)
    mol_3d = Chem.MolToMolBlock(mol)
    return mol_3d  


mol_3d = ligand_emb(mol)
print(mol_3d)

import pybel

for mol in pybel.readfile('sdf', 'many_molecules.sdf'):
  for atom in mol:
    coords = atom.coords
    for neighbor in pybel.ob.OBAtomAtomIter(atom.OBAtom):
      neighbor_coords = pybel.atom(neighbor).coords