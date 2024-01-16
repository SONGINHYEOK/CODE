import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from chembl_structure_pipeline import standardizer
from multiprocessing import Pool

def standardize_and_canonicalize_smiles(smiles):
    standardizer_result = standardizer.run(smiles)
    standardized_smiles = standardizer_result.standardized_smiles
    
    mol = Chem.MolFromSmiles(standardized_smiles)
    canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
    
    return canonical_smiles if canonical_smiles != smiles else smiles, mol

def process_row(row):
    canonical_smile, canonical_mol = standardize_and_canonicalize_smiles(row['smile'])
    
    # Set additional information as SD data
    sd_data = {'Index': str(row.name), 'AnotherColumnValue': str(row['another_column'])}
    for key, value in sd_data.items():
        canonical_mol.SetProp(key, value)
    
    return canonical_smile, canonical_mol

def process_csv(input_file, output_sdf, num_processes=4):
    # Read CSV
    df = pd.read_csv(input_file)

    # Use multiprocessing to process SMILES in parallel
    with Pool(num_processes) as pool:
        results = pool.map(process_row, df.itertuples(index=False))

    # Save all compounds in a single SDF file
    writer = Chem.SDWriter(output_sdf)
    for smile, mol in results:
        writer.write(mol)
    writer.close()

    # Update the 'smile' column with canonical SMILES
    df['smile'] = [smile for smile, _ in results]

    # Save the updated dataframe to a new CSV file
    df.to_csv("output_file.csv", index=False)

if __name__ == "__main__":
    input_csv = "your_input_file.csv"
    output_sdf = "output_file.sdf"
    num_processes = 4  # Adjust the number of processes as needed
    process_csv(input_csv, output_sdf, num_processes)
