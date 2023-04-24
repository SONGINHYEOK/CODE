import pandas as pd
import numpy as np
df_list = []
chunksize = 998
i = 0
df = pd.read_csv('/Users/song-inhyeok/Documents/BIChem/chembl_start_data_15000.csv')
print(len(df))
list_test = list(range(1,2326992000))


for chunk in np.array_split(df, len(df) // chunksize):
    chunk.to_csv(f"/Users/song-inhyeok/Desktop/admet_data/chembl_data/test/chembl_{i}.csv")
    #PandasTools.AddMoleculeColumnToFrame(chunk, smilesCol='smiles')
    #PandasTools.WriteSDF(chunk, f'/Users/song-inhyeok/Desktop/admet_data/chembl_data/test/chembl_{i}.sdf', molColName='ROMol', idName='id', properties=None, allNumeric=False)
    i = i + 1
    print(i)
