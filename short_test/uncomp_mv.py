import os


path = os.getcwd()

#폴더 내 tar.gz file list 만들기

gzfile_list = [path + "/conf_"+ str(i)+".tar.gz" for i ]