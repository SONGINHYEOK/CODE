import pandas as pd
from rdkit import Chem
len_list = []

df = pd.read_csv("/Users/song-inhyeok/CODING/short_test/not_cal_scaffold.csv")
smi_list = df['smiles'].tolist()

n = 0
for smi in smi_list:
    n +=1
    m = Chem.MolFromSmiles(smi)
    len_list.append(m.GetAtoms())

    if n == 10:
        break
print(len_list)