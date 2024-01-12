import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from multiprocessing import Pool, freeze_support

def process_sequence(sequence):
    mol = Chem.MolFromSequence(sequence)
    if mol is not None:
        AllChem.Compute2DCoords(mol)
        AllChem.EmbedMolecule(mol)
        return mol
    else:
        print(f"Failed to generate molecule from sequence: {sequence}")
        return None

def csv_to_sdf(csv_file, output_sdf, peptide_column, num_processes=4):
    # Read the CSV file and extract the peptide sequences
    df = pd.read_csv(csv_file)
    
    # Check if the specified peptide column exists
    if peptide_column not in df.columns:
        print(f"Error: '{peptide_column}' column not found in the CSV file.")
        return
    
    peptide_sequences = df[peptide_column]

    # Use multiprocessing to process peptide sequences in parallel
    with Pool(num_processes) as pool:
        mols = pool.map(process_sequence, peptide_sequences)

    # Filter out None values (failed sequences) and write valid molecules to an SDF file
    valid_mols = [mol for mol in mols if mol is not None]

    if valid_mols:
        # Write the molecules to an SDF file
        writer = Chem.SDWriter(output_sdf)
        for mol in valid_mols:
            writer.write(mol)
        writer.close()
        print(f"SDF file '{output_sdf}' created successfully.")
    else:
        print("No valid molecules to write.")

if __name__ == '__main__':
    freeze_support()  # Required for Windows support
    # Example usage:
    csv_file_path = './sample.csv'
    output_sdf_path = './output_conv_pep_molecules.sdf'
    peptide_column_name = 'Peptide Sequence'
    num_processes = 4  # Adjust the number of processes as needed
    csv_to_sdf(csv_file_path, output_sdf_path, peptide_column_name, num_processes)
