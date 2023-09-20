import multiprocessing
from rdkit import Chem
from rdkit.Chem import AllChem
from chembl_structure_pipeline import standardizer


# Function to generate conformers for a molecule
def generate_conformers(smiles, num_conformers, lock):
    mol = Chem.MolFromSmiles(smiles)
    s = Chem.MolToMolBlock(mol)
    std_molblock = standardizer.standardize_molblock(s)
    std_mol = Chem.MolFromMolBlock(std_molblock, removeHs=False)
    if std_mol is not None:
        AllChem.EmbedMultipleConfs(std_mol, numConfs=num_conformers)
        for conf_id in range(std_mol.GetNumConformers()):
            AllChem.MMFFOptimizeMolecule(std_mol, confId=conf_id)
        
        with lock:
            writer = Chem.SDWriter("/Users/song-inhyeog/CODEING/CODE/SMC_task/small_molecule_DB/molecule_perparation/confgen/conformers.sdf")
            
            for conf in std_mol.GetConformers():
                std_mol.SetProp('ID', f'{smiles}_{conf.GetId()}')
                writer.write(std_mol, confId=conf.GetId())
            writer.close()

def main():
    num_conformers = 10  # Number of conformers to generate
    num_processes = 2 #multiprocessing.cpu_count()  # Number of CPU cores

    # Create a manager and a shared lock
    with multiprocessing.Manager() as manager:
        lock = manager.Lock()
        
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
        
        
        # Split the conformer generation tasks across processes
        pool = multiprocessing.Pool(processes=num_processes)
        tasks = [(smiles, num_conformers, lock) for smiles in smiles_list]
        pool.starmap(generate_conformers, tasks)
        pool.close()
        pool.join()
        
        print(f"Conformers generated successfully for {len(smiles_list)} molecules.")

if __name__ == "__main__":
    main()
