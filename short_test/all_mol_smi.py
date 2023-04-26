import multiprocessing
import time
import os
import math
from multiprocessing import JoinableQueue
import shutil
import pandas as pd

def multimv(folder_list):
    start = time.time()
    #try:
    df = pd.read_csv(folder_list, header=None)
    smiles=df[0].tolist()
    ids = df[1].tolist()
    for smile,id in zip(smiles, ids):
        with open("./all_mol.smi", 'a') as f:
            f.write(str(smile) + ' ' + str(id) + "\n")
            f.close()
         
    

def worker(q):
    for item in iter(q.get, None):
        multimv(item)
        q.task_done()
    q.task_done()

def main():
    num_procs = 5
    path = "/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/all_chem_data"
    file_list = [path +'/all_mol_'+str(i)+'.csv' for i in range(1, 938)]
    #folder_list = [ path + '/conf_' + str(i) for i in range(1882, 2159)]
    q = multiprocessing.JoinableQueue()
    procs = []
    for i in range(num_procs):
        procs.append(multiprocessing.Process(target=worker, args=(q,) ))
        #proc = multiprocessing.Process(target=worker, args=(q,))
        #procs.append(proc)
        time.sleep(0.1)
        procs[-1].daemon = True
        procs[-1].start()
    for folder in file_list:
        q.put(folder)
    q.join()
    for p in procs:
        q.put(None)
    q.join()
    for p in procs:
        p.join()
if __name__ == '__main__':
    main()