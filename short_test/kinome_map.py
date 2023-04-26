import pandas as pd

df = pd.read_csv("/Users/song-inhyeok/CODING/short_test/kinome.txt", header=None, sep='\t')


node_list = []

for index, value in enumerate(df[0].tolist()):
     
    result = divmod(index, 2)
    
    if result[1] == 1:
    
        node_list.append(value)

name_list = []

for i in node_list:
    name_list.append(i.split(" ")[1].split("_")[2].replace('"',''))


result1 = set(name_list)
print(f"set(arr)      : {result1}")

result2 = list(result1)  # list(set(arr))
print(f"list(set(arr) : {result2}")

print(len(name_list))
print(len(result2))


df2 = pd.DataFrame()
df2['id']= result2

df2.to_csv("/Users/song-inhyeok/Desktop/kinomemap_labellist.csv")
