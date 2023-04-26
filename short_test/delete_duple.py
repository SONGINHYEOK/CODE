import pandas as pd

df = pd.read_csv("/Users/song-inhyeok/CODING/short_test/sch_not_ligand_list.text", sep='\t', header=None)
pdb_list=df[0].tolist()

result1 = set(pdb_list)
result2 = list(result1)  # list(set(arr))

print(len(pdb_list))
print(len(result1))
print(len(result2))


for i in result2:
    with open("./Sch_not_ligand_list.txt", 'a') as f:
        f.write(i + '\n')
        f.close()