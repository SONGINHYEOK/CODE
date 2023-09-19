import multiprocessing
from rdkit import Chem
from rdkit.Chem import AllChem

# Function to generate conformers for a molecule
def generate_conformers(smiles, num_conformers, lock):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers)
        for conf_id in range(mol.GetNumConformers()):
            AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)
        
        with lock:
            writer = Chem.SDWriter("conformers.sdf")
            
            for conf in mol.GetConformers():
                mol.SetProp('ID', f'{smiles}_{conf.GetId()}')
                writer.write(mol, confId=conf.GetId())
            writer.close()

def main():
    num_conformers = 10  # Number of conformers to generate
    num_processes = 4 #multiprocessing.cpu_count()  # Number of CPU cores

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
