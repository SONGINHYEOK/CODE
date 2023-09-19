from rdkit import Chem
from rdkit.Chem import AllChem

# Define your molecule as a SMILES string or load it from a file
smiles = "CCO"  # Example SMILES for ethanol

# Convert the SMILES string to an RDKit molecule object
mol = Chem.MolFromSmiles(smiles)

if mol is not None:
    # Generate conformers
    num_conformers = 10  # Number of conformers to generate
    AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers)

    # You can also optimize the generated conformers
    for conf_id in range(mol.GetNumConformers()):
        AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)

    # Save the conformers to an SDF file
    writer = Chem.SDWriter("conformers.sdf")
    for conf in mol.GetConformers():
        writer.write(mol, confId=conf.GetId())

    writer.close()

    # Alternatively, you can access and manipulate individual conformers
    conformer_0 = mol.GetConformer(0)
    # Do something with the conformer_0

else:
    print("Failed to generate molecule from SMILES.")
