import tarfile
import time
import os
import math
import shutil

## 1 file extract & each file move
start = time.time()

# make file path 
file_list = ['/home/bichem_storage/conformer/conf_' + str(i)+".tar.gz" for i in range(4403,4582)]
#file_path = '/Users/johyungsik/Code/PDB/conf_5000.tar.gz'

# extract

target_path = '/home/bichem_storage/conformer/'
for file_path in file_list:
    tar = tarfile.open(file_path, 'r:gz')

    for item in tar:
        #count += 1
        one_time = time.time()
        tar.extract(item, target_path)

        file_num = item.name.split('.')[0].split('/')[-1]
        file_name = item.name.split('/')[-1]

        old_folder_path = target_path + item.name
        new_folder_path = target_path + 'conformer_' + str(math.floor(int(file_num)/500000)) + '/'

        if not os.path.isdir(new_folder_path):
            os.mkdir(new_folder_path)
            
        new_file_path = new_folder_path + file_name
        shutil.move(old_folder_path, new_file_path)

        #one_end_time = time.time() - one_time
        #print('1 file extract & move time : ', one_end_time, ', move_path : ', new_file_path)
    tar.close()

end = time.time() - start
print('end time', end)