






import os
import tarfile
import math


root_path = os.getcwd()
conf_path = [root_path + "/conf_" + str(i)+".tar.gz" for i in range(4620,4621)]


for file in conf_path:
    mid_path = file.split(".")[0] + "/"
    tar = tarfile.open(file)
    tar_list = tar.getnames()
    #if get file list
    if len(tar_list) > 10000000:
        os.system("tar -zxf"+ " "+ str(file))
        for tar_file in tar_list:
        pre_num=tar_file.split('.')[0].split("/")[1]
        num = int(math.floor(int(pre_num)/500000))
        print( "file mv")
        os.system("mv"+" "+ root_path+ "/" + tar_file+" "+"/home/storage/conformer/"+"conformer_" + str(num))
    
    else:
        os.system("tar -zxf"+ " "+ str(file))




import os

root_path = os.getcwd()
conf_path = [root_path + "/conf_" + str(i)+".tar.gz" for i in range(4621,4625)]

for file in conf_path:
    os.system("tar -zxf"+ " "+ str(file))