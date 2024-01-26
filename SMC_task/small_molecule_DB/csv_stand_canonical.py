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
    sd_data = {'Index': str(row['Drug_ID']), 'target': str(row['Y'])}
    for key, value in sd_data.items():
        canonical_mol.SetProp(key, value)
    
    return {'Drug_ID': row['Drug_ID'], 'Y': row['Y'], 'canonical_smile': canonical_smile, 'canonical_mol': canonical_mol}

def process_csv(input_file, output_sdf, output_csv, num_processes=4):
    # Read CSV
    df = pd.read_csv(input_file)

    # Use multiprocessing to process SMILES in parallel
    with Pool(num_processes) as pool:
        results = pool.map(process_row, df.itertuples(index=False))

    # Save all compounds in a single SDF file
    writer = Chem.SDWriter(output_sdf)
    for result in results:
        canonical_mol = result['canonical_mol']
        writer.write(canonical_mol)
    writer.close()

    # Update the 'smile' column with canonical SMILES
    df['smile'] = [result['canonical_smile'] for result in results]

    # Save the updated dataframe to a new CSV file
    df.to_csv(output_csv, index=False, columns=['Drug_ID', 'Y', 'smile'])

if __name__ == "__main__":
    input_csv = "/Users/song-inhyeog/Downloads/caco2_total.csv"
    output_sdf = "output_std.sdf"
    output_csv = "output_std_file.csv"
    num_processes = 4  # Adjust the number of processes as needed
    process_csv(input_csv, output_sdf, output_csv, num_processes)
