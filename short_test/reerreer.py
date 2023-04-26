import multiprocessing
import time
import os
import math
from multiprocessing import JoinableQueue
import shutil
#동작 실행 함수

def multimv(folder_list):
    for filename in os.listdir(folder_list):
        file_path =filename
        pre_num=filename.split('.')[0]
        
        old_folder_path = folder + file_path
        new_folder_path =  '/home/bichem_storage/conformer/'+ 'conformer_' + str(math.floor(int(pre_num)/500000)) + '/'

        if not os.path.isdir(new_folder_path):
            os.mkdir(new_folder_path)
            
        new_file_path = new_folder_path + str(file_path)
        shutil.move(old_folder_path, new_file_path)
        

#Queue 전달 함수
def worker(q):
    for item in iter(q.get, None):
        multimv(item)
        q.task_done()
    q.task_done()

#메인 함수
def main():

    #동작 프로세스 개수
    num_procs = 4

    #큐 데이터

    path = os.getcwd()
    folder_list = [ path + '/conf_' + str(i) for i in range(240, 1801)]

    q = multiprocessing.JoinableQueue()
    procs = []
    
    for i in range(num_procs):
        procs.append(multiprocessing.Process(target=worker, args=(q,) ))
        #proc = multiprocessing.Process(target=worker, args=(q,))
        #procs.append(proc)
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

    #q.__init__(ctx=None)
    #q.__init__()
    print ("Finished everything....")
    print ("num active children:", multiprocessing.active_children())


if __name__ == '__main__':
    main()