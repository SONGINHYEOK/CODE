import gzip

file_path = '/Users/song-inhyeog/CODEING/CODE/SMC_task/small_molecule_DB/molecule_perparation/confgen/H29P900.smi.gz'

with gzip.open(file_path, 'rt') as gz_file:
    content = gz_file.read().split('\n')
    del content[-1]

smiles_list = []
for i in content:
    smi = i.split(" ")[0].split("\t")[0]
    #ori_id = i.split(" ")[0].split("\t")[1]
    smiles_list.append(smi)

#=====================================================================================================#

import os
import multiprocessing
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import SDWriter

# Define a function to generate conformers for a single molecule
def generate_conformers(molecule, num_conformers):
    mol = Chem.MolFromSmiles(molecule)
    if mol is None:
        return None
    AllChem.EmbedMultipleConfs(mol, numConformers=num_conformers)
    return mol

# Define a function to process a batch of molecules in parallel
def process_batch(batch, num_conformers, writer):
    for molecule in batch:
        conformer = generate_conformers(molecule, num_conformers)
        if conformer:
            writer.write(conformer)  # Write the conformer to the SDF file

if __name__ == '__main__':
    # Define your list of SMILES strings for the molecules you want to generate conformers for
    

    # Specify the number of conformers to generate for each molecule
    num_conformers_per_molecule = 20

    # Create an SDF writer to save the conformers
    sdf_filename = 'conformers.sdf'
    writer = SDWriter(sdf_filename)

    # Split the list of SMILES into batches for parallel processing
    num_processes = 4#multiprocessing.cpu_count()
    batch_size = len(smiles_list) // num_processes
    batches = [smiles_list[i:i+batch_size] for i in range(0, len(smiles_list), batch_size)]

    # Create a multiprocessing pool to generate conformers in parallel
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.starmap(process_batch, [(batch, num_conformers_per_molecule, writer) for batch in batches])

    # Close the SDF writer
    writer.close()

    print(f'Conformers saved to {sdf_filename}')
