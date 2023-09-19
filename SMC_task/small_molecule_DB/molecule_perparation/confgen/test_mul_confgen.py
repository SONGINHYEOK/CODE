import multiprocessing
from rdkit import Chem
from rdkit.Chem import AllChem

# Define your molecule as a SMILES string or load it from a file
smiles = "CCO"  # Example SMILES for ethanol

def generate_conformers(mol, num_conformers, lock):
    """Generate conformers for a molecule."""
    AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers)
    
    # Optimize the generated conformers
    for conf_id in range(mol.GetNumConformers()):
        AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)
    
    # Write conformers to a file with a lock to prevent race conditions
    with lock:
        writer = Chem.SDWriter("conformers.sdf")
        for conf in mol.GetConformers():
            writer.write(mol, confId=conf.GetId())
        writer.close()

def main():
    # Convert the SMILES string to an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)

    if mol is not None:
        num_conformers = 10  # Number of conformers to generate
        num_processes = multiprocessing.cpu_count()  # Number of CPU cores
        
        # Create a lock to protect file writing
        lock = multiprocessing.Lock()
        
        # Split the conformer generation tasks across processes
        pool = multiprocessing.Pool(processes=num_processes)
        tasks = [(mol, num_conformers, lock)] * num_processes
        pool.starmap(generate_conformers, tasks)
        pool.close()
        pool.join()
        
        print(f"Conformers generated successfully ({num_conformers * num_processes} conformers).")

    else:
        print("Failed to generate molecule from SMILES.")

if __name__ == "__main__":
    main()