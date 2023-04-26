import os
import math
import time

     
path = os.getcwd()
folder_list = [ path + '/conf_' + str(i) for i in range(1)]


for folder in folder_list:
    for i in range(len(os.listdir(folder))):
        one_start = time.time()
        file_path =os.listdir(folder)[i]
        pre_num=os.listdir(folder)[i].split('.')[0]
        num = int(math.floor(int(pre_num)/500000))
        os.system("mv"+" "+folder+ "/" +str(file_path) +" "+"/home/bichem_storage/conformer/"+"conformer_" + str(num))
        print(file_path , ' conformer_' , num)
        print("time :", time.time() - one_start)    
        
        if i == 9:
            break

