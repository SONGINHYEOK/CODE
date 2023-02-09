from rdkit import Chem, RDConfig
from rdkit.Chem import MolStandardize, Draw, AllChem, rdMolAlign
from rdkit.Chem import PandasTools
import pandas as pd
name_list = []

file_name="/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/DOWNLOAD--6WXs2QhCLCIkPW3JobW_iub6DxpLqtZK4IWp5Xzd5s=.sdf"
df=PandasTools.LoadSDF(file_name)

#print(df)
"""
sdf = Chem.SDMolSupplier(file_name)
for i in sdf:
    #score_list.append(i.GetProp('r_i_docking_score'))
    print(i.GetProp('chembl_id'))
"""
id_list = []
ori_list=df['chembl_id'].tolist()
for i in ori_list: 
    ori_id=i[6:]
    id_list.append(ori_id)
    
df2 = pd.DataFrame()  
df2['id']=id_list
df2.to_csv("./appoped_id.csv", index=False)
