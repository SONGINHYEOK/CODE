from os import remove
from matplotlib.cbook import sanitize_sequence
from rdkit import Chem, RDConfig
from rdkit.Chem import Draw, AllChem, rdMolAlign, MolStandardize
import time


def _InitialiseNeutralisationReactions():
    patts= (
        # Imidazoles
        ('[n+;H]','n'),
        # Amines
        ('[N+;!H0]','N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]','O'),
        # Thiols
        ('[S-;X1]','S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]','N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]','N'),
        # Tetrazoles
        ('[n-]','[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]','S'),
        # Amides
        ('[$([N-]C=O)]','N'),
        )
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in patts]

_reactions=None
def NeutraliseCharges(smiles, reactions=None):
    global _reactions
    if reactions is None:
        if _reactions is None:
            _reactions=_InitialiseNeutralisationReactions()
        reactions=_reactions
    mol = Chem.MolFromSmiles(smiles)
    replaced = False
    for i,(reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant):
            replaced = True
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]
    if replaced:
        return (Chem.MolToSmiles(mol,True), True)
    else:
        return (smiles, False)



#mols = Chem.SDMolSupplier('/Users/song-inhyeok/Documents/data/all-sdf.sdf')
mols = Chem.SDMolSupplier('/Users/song-inhyeok/Documents/data/all-sdf.sdf',sanitize = False)

from rdkit import Chem
from rdkit.Chem import AllChem



from openbabel import pybel

for mol in pybel.readfile('sdf','/Users/song-inhyeok/Documents/data/all-sdf.sdf'):
  for atom in mol:
    coords = atom.coords
    for neighbor in pybel.ob.OBAtomAtomIter(atom.OBAtom):
      neighbor_coords = pybel.atom(neighbor).coords





start = time.time()

def ligand_emb(mol):
    mol2 = Chem.AddHs(mol)
    AllChem.ConstrainedEmbed(mol2, mol)
    mol_3d = Chem.MolToMolBlock(mol)
    return mol_3d  

for mol in mols:
    try:
        ligand_3d = ligand_emb(mol)
    except:
        continue    
        
        
    except:
        with open("./error.txt", "a") as f:
            f.write(f"{}\n") 
            f.close()  
    try:
        (molSmiles, neutralised) = NeutraliseCharges(ligand_smi)
        standard_ligand_smi = MolStandardize.canonicalize_tautomer_smiles(molSmiles)
        print(standard_ligand_smi)
    except IndexError:
        print('{} index의 값을 가져올 수 없습니다.'.format(standard_ligand_smi))
    except AttributeError:
        print('{} attribute의 값을 가져올 수 없습니다.'.format(standard_ligand_smi))

    
end = time.time() - start
print('Running time;', end)
print("Job done")    
        
   
#!/bin/bash

for i in {3291..3721}
do
    tar -zxf conf_$i.tar.gz
done