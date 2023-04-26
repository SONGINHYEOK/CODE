import math
import time
import subprocess
import os


def mv_func(tar_list, mid_path):
    for tar_file in tar_list:
        pre_num=tar_file.split('.')[0].split("/")[1]
        num = int(math.floor(int(pre_num)/500000))
        os.system("mv"+" "+ mid_path  + tar_file+" "+"/home/storage/conformer/"+"conformer_" + str(num))
        
if __name__ == '__main__':
    mv_func()