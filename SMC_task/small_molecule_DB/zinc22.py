import multiprocessing
import time
import os
import math
from multiprocessing import JoinableQueue
import shutil
import pandas as pd

import urllib.request as urllib


def multimv(folder_list):
    count = 0
    start = time.time() 
    down_path = "/Users/song-inhyeog/CODEING/CODE/SMC_task/small_molecule_DB/zinc22_2d/original_file/"
    file_name=folder_list.split('/')[-1]
    
 
    final_path = down_path + file_name
    urllib.urlretrieve(folder_list, final_path)
    #try:
    #df = pd.read_csv('/Users/song-inhyeok/Documents/ZINC20_cr_project/ZINC20_smi/0/'+folder_list, header=None)
    #zid=df[0].tolist()
    #print(zid)
    #for i in zid:
    #    url = 'https://zinc20.docking.org/substances/' + i
    #   try:
    #        urllib.urlretrieve(url, '/Users/song-inhyeok/Documents/ZINC20_cr_project/ZINC_ALL/'+i)
    #    except:
    #        continue
    #    count +=1
    #except:
    #    file_list = []
    #    print("don't have folder : ", folder_list)
 
    
def worker(q):
    for item in iter(q.get, None):
        multimv(item)
        q.task_done()
    q.task_done()

def main():
    num_procs = 4
    
    f = open("/Users/song-inhyeog/CODEING/CODE/SMC_task/small_molecule_DB/zinc22_2d/download_link/ZINC22_12.txt", "r")
    
    lst = []
    while True:
        line = f.readline().strip()
        if not line:break
        lst.append(line)
    
    
    q = multiprocessing.JoinableQueue()
    procs = []
    for i in range(num_procs):
        procs.append(multiprocessing.Process(target=worker, args=(q,) ))
        #proc = multiprocessing.Process(target=worker, args=(q,))
        #procs.append(proc)
        time.sleep(0.1)
        procs[-1].daemon = True
        procs[-1].start()
    for folder in lst:
        q.put(folder)
    q.join()
    for p in procs:
        q.put(None)
    q.join()
    for p in procs:
        p.join()
if __name__ == '__main__':
    main()