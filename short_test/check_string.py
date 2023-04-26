import pandas as pd

id_list = []
smi_list = []


#ori_df = pd.read_csv("/Users/song-inhyeok/Desktop/admet_data/chembl_data/chembl_start_data.csv")


#ids = ori_df['id'].tolist()
#smiles=ori_df["smiles"].tolist()
with open("/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/sure_chembl.csv", "r") as f:
    line=f.readlines()
    for li in line:
        print(li)
        """
        ids=li.spilt(',')[0]
        smiles=li.splilt(",")[1]
for id, smile in zip(ids,smiles):
    if len(smile) > 415 :
        id_list.append(id)
        smi_list.append(smile)

df = pd.DataFrame()
df['id'] = id_list
df['smiles']= smi_list

df.to_csv("./not_cal_surechem.csv", index=False)
"""