import pandas as pd
import numpy as np
zid_list = []
chunksize = 10000
n = 0
#zinc_20 = ['ZINC'+str(i).zfill(12) for i in range(1,2326992001)]
#zinc_20 = ['ZINC'+str(i).zfill(12) for i in range(1,2326992001)]


for i in range(100000000, 200000000):
    zid = 'ZINC'+str(i).zfill(12)+'.smi'
    zid_list.append(zid)


df = pd.DataFrame(zid_list)
#df.to_csv('./total_0', index=False, header=False)
#zinc_20 = ['ZINC'+str(i).zfill(12) for i in range(1,23)]
#print(zinc_20[:10])
#list_chunked = list_chunk(zinc_20, 116349600) 

#print(list_chunked[0])
j = 1
for chunk in np.array_split(df, len(df) // chunksize):
    chunk.to_csv(f"/Users/song-inhyeok/Documents/ZINC20_cr_project/ZINC20_smi/{j}/ZINC_{j}_{n}.csv", index=False, header=False)
    n +=1
    print(n)