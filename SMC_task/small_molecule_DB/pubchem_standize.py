from rdkit import Chem
from rdkit.Chem import AllChem
from chembl_structure_pipeline import standardizer
import pandas as pd
def read_sdf_and_standardize(input_sdf, output_sdf):
    # Load the SDF file
    suppl = Chem.SDMolSupplier(input_sdf)

    # Initialize a list to store standardized molecules
    standardized_molecules = []

    # Process each molecule in the SDF file
    for mol in suppl:     
        if mol is not None:
            try:
                as_mol = AllChem.AssignStereochemistry(mol)
                ach_mol = AllChem.AssignAtomChiralTagsFromStructure(as_mol)
                standardized_molecules.append(ach_mol)
            except:
                standardized_molecules.append(mol)
        
    smi_list = []
    for mol in standardized_molecules:
        smi = Chem.MolToSmiles(mol)
        smi_list.append(smi)
        
    df = pd.DataFrame()
    df['smile'] = smi_list
    df.to_csv('pubchem_std.csv')
    # Write the standardized molecules to a new SDF file
    #w = Chem.SDWriter(output_sdf)
    #for mol in standardized_molecules:
    #    w.write(mol)
    #w.close()

if __name__ == "__main__":
    # Provide the input and output SDF file paths
    input_sdf_file = "/home/song/Documents/CODING/data/test/drugbank_approved.sdf"
    output_sdf_file = "output_pubchem_standardized.sdf"

    # Call the function to read SDF, convert to SMILES, and standardize
    read_sdf_and_standardize(input_sdf_file, output_sdf_file)
