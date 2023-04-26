import pandas as pd

df = pd.read_csv("/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/chembl.csv", header=None)

id_list = df[0].tolist()

ori_df = pd.read_csv("/Users/song-inhyeok/Downloads/chembl_20230116.csv", sep=';')

print(df.info())



# mask에다가 어떤 칼럼명에 어떤 값을 갖고 있는 애들을 선택하게 한다.

mask = ori_df["ChEMBL ID"].isin(id_list)
 

# 저장하려면 객체에 다시 저장

 

finanl_df = ori_df[~mask]

finanl_df = finanl_df.drop_duplicates(["ChEMBL ID"]) 
print(len(finanl_df["ChEMBL ID"].tolist()))


finanl_df.to_csv("/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/not_add_chembl.csv",index=False)
 
 

