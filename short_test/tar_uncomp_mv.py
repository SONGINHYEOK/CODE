import os
import tarfile
import math


root_path = os.getcwd()
conf_path = [root_path + "/conf_" + str(i)+".tar.gz" for i in range(4620,4621)]


for idx,file in enumerate(conf_path):
    tar = tarfile.open(file)
    tar_list = tar.getnames()
    #if get file list
    if len(tar_list) > 10000000:
        os.system("tar -zxf"+ " "+ str(file))
        for tar_file in tar_list:            
            pre_num=tar_file.split('.')[0]    
            num = int(math.floor(int(pre_num)/500000))
            file_path = root_path + "/" + "con_f"+str(num) + tar_file
            os.system("cp"+" "+ str(file_path) +" "+"/home/bichem_storage/conformer/"+"conformer_" + str(num))
                
    else:
        os.system("tar -zxf"+ " "+ str(file))
        
        
        
        
        


