import os
import tarfile
import math


path = "/Users/song-inhyeok/Documents/data/tar_test/tar_total2"
conf_path = [ path + '/test_' + str(i)+".tar.gz" for i in range(10)]

#root_path = os.getcwd()
#conf_path = [root_path + "/conf_" + str(i)+".tar.gz" for i in range(3764,3765)]


for idx,file in enumerate(conf_path):
    print(file)
    tar = tarfile.open(file)
    print(tar)
    tar_list = tar.getnames()
    tar.extractall()
    tar.close()
    #if get file list
    if len(tar_list) > 4000:
        for tar_file in tar_list:            
            pre_num=tar_file.split('.')[0]    
            num = int(math.floor(int(pre_num)/500000))
            #file_path = root_path + "/" + "con_f"+str(num) + tar_file
            os.system("cp"+" "+path+ "/" "test_"+str(idx) +" "+"/Users/song-inhyeok/Documents/data/tar_test/test_total2/" )
    else:
        os.system("tar -zxvf"+ " "+ str(file))
        
        
        
        
        


