from rdkit import Chem
from rdkit.Chem import rdRGroupDecomposition
from rdkit.Chem import PandasTools
import pandas as pd
import numpy as np


chunksize = 1000
i = 0
df = pd.read_csv("/Users/song-inhyeok/Desktop/admet_data/chembl_data/chembl_start_data.csv")

for chunk in np.array_split(df, len(df) // chunksize):
    PandasTools.AddMoleculeColumnToFrame(chunk, smilesCol='smiles')
    PandasTools.WriteSDF(chunk, f'/Users/song-inhyeok/Desktop/admet_data/chembl_data/sdf_file/chembl_{i}.sdf', molColName='ROMol', idName='id', properties=None, allNumeric=False)
    chunk.to_csv(f"/Users/song-inhyeok/Desktop/admet_data/chembl_data/csv_file/chembl_{i}.csv")
    i = i + 1
