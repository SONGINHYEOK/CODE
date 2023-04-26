import pandas as pd

df = pd.read_csv("/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/not_scaf_surechem.csv", header=None)
du_list = df[0].tolist()
print(len(du_list))
result1 = set(du_list)


result2 = list(result1)  # list(set(arr))



print(len(result2))


df2 = pd.DataFrame()
df2['id'] = result2
df2.to_csv("/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/not_scaf_surechem_del_duple.csv", index=False)