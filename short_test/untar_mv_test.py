import tarfile
import asyncio
import time
import os
import math

root_path = "/Users/song-inhyeok/Documents/data/tar_test/tar_total2"
conf_path = [root_path + "/test_" + str(i)+".tar.gz" for i in range(10)]


async def filemv_async(tar_list):
    for tar_file in tar_list:
        print(tar_file)    
        pre_num=tar_file.split('.')[0]
        num = int(math.floor(int(pre_num)/500000))
        print(num)
        os.system("mv"+" "+ root_path  + tar_file+" "+"/Users/song-inhyeok/Documents/data/tar_test/test_total2/")
        
async def main():
    for file in conf_path:
        print(file)
        tar = tarfile.open(file)
        tar_list =tar.getnames()
        tar.extractall(root_path)
        asyncio.sleep(5)
        await filemv_async(tar_list)


if __name__ == '__main__':
    asyncio.run(main())
    #asyncio.get_event_loop().run_until_complete(main())