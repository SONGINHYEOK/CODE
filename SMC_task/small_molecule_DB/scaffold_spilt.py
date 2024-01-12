from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import SDMolSupplier
from sklearn.model_selection import train_test_split
import pandas as pd

def generate_scaffold(molecule):
    scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
    return Chem.MolToSmiles(scaffold)

def split_dataset(sdf_file, test_size=0.2, valid_size=0.2, random_seed=42):
    # Read molecules from SDF file
    suppl = SDMolSupplier(sdf_file)

    # Generate scaffolds for each molecule
    scaffolds = [generate_scaffold(mol) for mol in suppl if mol is not None]

    # Split into train, test, and valid sets based on scaffolds
    scaffold_train, scaffold_test = train_test_split(scaffolds, test_size=test_size, random_state=random_seed)
    scaffold_train, scaffold_valid = train_test_split(scaffold_train, test_size=valid_size, random_state=random_seed)

    # Collect indices of molecules for each set
    train_indices = [i for i, mol in enumerate(suppl) if mol is not None and generate_scaffold(mol) in scaffold_train]
    test_indices = [i for i, mol in enumerate(suppl) if mol is not None and generate_scaffold(mol) in scaffold_test]
    valid_indices = [i for i, mol in enumerate(suppl) if mol is not None and generate_scaffold(mol) in scaffold_valid]

    return train_indices, test_indices, valid_indices

# Example usage
sdf_file_path = '/Users/song-inhyeog/CODEING/CODE/SMC_task/peptide/output_conv_pep_molecules.sdf'
train_indices, test_indices, valid_indices = split_dataset(sdf_file_path)

# Now, you can access the molecules in each set using the indices
train_molecules = [mol for i, mol in enumerate(SDMolSupplier(sdf_file_path)) if i in train_indices]
test_molecules = [mol for i, mol in enumerate(SDMolSupplier(sdf_file_path)) if i in test_indices]
valid_molecules = [mol for i, mol in enumerate(SDMolSupplier(sdf_file_path)) if i in valid_indices]



# Extract SMILES for each valid molecule
train_smiles = [{'ID': i, 'SMILES': Chem.MolToSmiles(mol)} for i, mol in enumerate(train_molecules) if i in valid_indices]

# Create a DataFrame with ID and SMILES columns for valid molecules
train_df = pd.DataFrame(train_smiles)

# Save the DataFrame to a CSV file
train_df.to_csv('train_set.csv', index=False)


# Extract SMILES for each valid molecule
test_smiles = [{'ID': i, 'SMILES': Chem.MolToSmiles(mol)} for i, mol in enumerate(test_molecules) if i in valid_indices]

# Create a DataFrame with ID and SMILES columns for valid molecules
test_df = pd.DataFrame(test_smiles)

# Save the DataFrame to a CSV file
test_df.to_csv('test_set.csv', index=False)


# Extract SMILES for each valid molecule
valid_smiles = [{'ID': i, 'SMILES': Chem.MolToSmiles(mol)} for i, mol in enumerate(valid_molecules) if i in valid_indices]

# Create a DataFrame with ID and SMILES columns for valid molecules
valid_df = pd.DataFrame(valid_smiles)

# Save the DataFrame to a CSV file
valid_df.to_csv('valid_set.csv', index=False)