import os
from Bio.PDB import PDBParser, PDBIO, Select
import Bio.PDB
import time

def is_het(residue):
    res = residue.id[0]
    return res != " " and res != "W"

class NonHetSelect(Select):
    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0

class ResidueSelect(Select):
    def __init__(self, chain, residue):
        self.chain = chain
        self.residue = residue

    def accept_chain(self, chain):
        return chain.id == self.chain.id

    def accept_residue(self, residue):
        """ Recognition of heteroatoms - Remove water molecules """
        return residue == self.residue and is_het(residue)


def extract_ligands(path):
    """ Extraction of the heteroatoms of .pdb files """
    i=1
    for pfb_file in os.listdir(path):
        if pfb_file.endswith('.pdb') and not pfb_file.startswith("lig_"):
            pdb_code = pfb_file[:-4]
            pdb = PDBParser().get_structure(pdb_code, path + pfb_file)
            io = PDBIO()
            io.set_structure(pdb)
            for model in pdb:
                model_atoms = Bio.PDB.Selection.unfold_entities(model, 'A')
                io.save(path+f"/{pdb_code}_rec.pdb", NonHetSelect())
                for chain in model:
                    print(chain)
                    for residue in chain:
                        if not is_het(residue):
                            continue
                        #print(residue.id)
                        print(f"saving {chain.id[0]} {residue}")
                        name = residue.id[0].split("_")[1]
                        
                        io.save(path+'/ligand/'+f"{pdb_code}_{chain.id[0]}_{name}_{residue.id[1]}_{i}.pdb", ResidueSelect(chain, residue))
                        i +=1
                        #io.save(path+'/ligand/'+f"{pdb_code}_{residue.id[0]}_{i}.pdb", ResidueSelect(chain, residue))
                                

# Main
path = '/Users/song-inhyeok/CODING/short_test/PDB/'

extract_ligands(path)




root="/Users/song-inhyeok/CODING/short_test/PDB/ligand/"

from openbabel import openbabel
import os, sys, copy
from openbabel import pybel

ligand_list = os.listdir("/Users/song-inhyeok/CODING/short_test/PDB/ligand")

obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("pdb", "sdf")

mol = openbabel.OBMol()

for ligand in ligand_list:
    name  = ligand.split(".")[0]
     
    
    obConversion.ReadFile(mol, root+ligand) # Open Babel will uncompress automatically
    
    mol.AddHydrogens()

    obConversion.WriteFile(mol,root+name+".sdf")




