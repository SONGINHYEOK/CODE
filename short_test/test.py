i = 0
with open ("/Users/song-inhyeok/Documents/ZINC20_cr_project/final_filtering_ZINC20.smi") as f:
    lines=f.readlines()
    for line in lines:
        print(line)
        i +=1
        if i==2 :
            break