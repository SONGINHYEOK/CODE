from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
import pandas as pd
from rdkit import Chem
from multiprocessing import Pool
from functools import partial

# Function to calculate descriptors for a single molecule
def calculate_descriptors(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        return calc.CalcDescriptors(mol)
    else:
        return [None] * len(desc_list)

# Function to create molecule from SMILES string
def mol_from_smiles(smi):
    return Chem.MolFromSmiles(smi)

# Number of processes to use
num_processes = 8  # Adjust this based on your system's capabilities

# Read the CSV file
start_df = pd.read_csv("/home/song/Documents/CODING/pubchem_std.csv")
smi_list = start_df['smile'].to_list()

# Create a pool of processes
with Pool(processes=num_processes) as pool:
    mols = pool.map(mol_from_smiles, smi_list)

# Create a descriptor calculator
desc_list = [n[0] for n in Descriptors._descList]
calc = MoleculeDescriptors.MolecularDescriptorCalculator(desc_list)

# Calculate descriptors using multiple processes
with Pool(processes=num_processes) as pool:
    rdkit_desc = pool.map(calculate_descriptors, smi_list)

# Create DataFrames
final_df = pd.DataFrame(rdkit_desc, columns=desc_list)
final_df['ID'] = range(len(smi_list))
final_df.to_csv("molecular_all_descriptor.csv", index=False)

info_df = pd.DataFrame()
info_df['ID'] = range(len(smi_list))
info_df['smile'] = smi_list
info_df.to_csv("ori_info.csv", index=False)