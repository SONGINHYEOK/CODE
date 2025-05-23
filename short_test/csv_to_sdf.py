from rdkit import Chem
from rdkit.Chem import rdRGroupDecomposition
from rdkit.Chem import PandasTools
import pandas as pd
import numpy as np


chunksize = 1000
i = 0
df = pd.read_csv("/Users/song-inhyeok/Desktop/admet_data/chembl_data/chembl_start_data.csv")

PandasTools.AddMoleculeColumnToFrame(df, smilesCol='smiles')
PandasTools.WriteSDF(df, './chembl_all.sdf', molColName='ROMol', idName='id', properties=None, allNumeric=False)
