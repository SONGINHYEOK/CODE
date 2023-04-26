import os
import math
import time

     
path = "/Users/song-inhyeok/Documents/data/tar_test"
folder_list = [ path + '/test_' + str(i) for i in range(10)]


start = time.time()

for folder in folder_list:
        for idx,filename in enumerate(os.listdir(folder)):
            one_start = time.time()
            file_path =filename
            pre_num=filename.split('.')[0]
            num = int(math.floor(int(pre_num)/500000))
            os.system("cp"+" "+folder+ "/" +str(file_path) +" "+"/Users/song-inhyeok/Documents/data/tar_test/test_total/" )
            print(file_path , ' conformer_' , num)
            print("time :", time.time() - one_start)
            
            
print("time :", time.time() - start)
            