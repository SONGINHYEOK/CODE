import pandas as pd
import urllib.request

save_path = '/Users/song-inhyeok/Documents/data/protein_img/'

df = pd.read_csv("/Users/song-inhyeok/Downloads/rcsb_pdb_ids_d86e7fdf4d640e40539a553fb8974089_100001-125000.txt", sep=',')

count = 100000
for i,name in enumerate(df.columns.tolist()):
    count +=1
    full_name = name.lower()
    mid_name = name.lower()[1:3]
    base_url = "https://cdn.rcsb.org/images/structures/"
    img_url = f'{base_url}/{mid_name}/{full_name}/{full_name}_assembly-1.jpeg'
    img_name = f'{name}.jpeg'
    urllib.request.urlretrieve(img_url, save_path+img_name)
    print(count)