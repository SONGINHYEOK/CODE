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
    #try:
    df = pd.read_csv('/Users/song-inhyeok/Documents/ZINC20_cr_project/ZINC20_smi/19/'+folder_list, header=None)
    zid=df[0].tolist()
    #print(zid)
    for i in zid:
        url = 'https://zinc20.docking.org/substances/' + i
        try:
            urllib.urlretrieve(url, '/Users/song-inhyeok/Documents/ZINC20_cr_project/ZINC_ALL/'+i)
        except:
            continue
        count +=1
    #except:
    #    file_list = []
    #    print("don't have folder : ", folder_list)
    end = time.time() - start
    print('end time : ', end, 'total count :', count)
  
    
def worker(q):
    for item in iter(q.get, None):
        multimv(item)
        q.task_done()
    q.task_done()

def main():
    num_procs = 10
    
    file_path = '/Users/song-inhyeok/Documents/ZINC20_cr_project/ZINC20_smi/19/'
    folder_list=os.listdir(file_path)
    q = multiprocessing.JoinableQueue()
    procs = []
    for i in range(num_procs):
        procs.append(multiprocessing.Process(target=worker, args=(q,) ))
        #proc = multiprocessing.Process(target=worker, args=(q,))
        #procs.append(proc)
        time.sleep(0.1)
        procs[-1].daemon = True
        procs[-1].start()
    for folder in folder_list:
        q.put(folder)
    q.join()
    for p in procs:
        q.put(None)
    q.join()
    for p in procs:
        p.join()
if __name__ == '__main__':
    main()