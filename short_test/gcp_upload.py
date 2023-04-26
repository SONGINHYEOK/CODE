from datetime import datetime, timedelta
#from bs4 import BeautifulSoup
from selenium import webdriver
from selenium.webdriver.support.ui import Select
from time import sleep
from os import listdir
from os.path import isfile, join
import os
import subprocess

from google.cloud import storage
os.environ["GLOUD_PROJECT"] = "nettargets-chemdb-v1"
#client = storage.Client()z
client = storage.Client.from_service_account_json('/Users/song-inhyeok/Documents/data/nettargets-chemdb-v1-beb86f192c43.json')
bucket = client.get_bucket('nettargets-discoveryweb-protein-bucket')
# blob = bucket.list_blobs('PDB/')
blob = storage.Blob('PDB/', bucket)


unzip_files = [f for f in listdir('/Users/song-inhyeok/Documents/coding/PDB/PDB/') if isfile(join('/Users/song-inhyeok/Documents/coding/PDB/PDB', f))]
unzip_files = [x for x in unzip_files if x.find("pdb") != -1]

print(unzip_files)






for i in unzip_files:
    try:
        file_name_rm = "/Users/song-inhyeok/Documents/coding/PDB/PDB/" + i
        blob.upload_from_filename(file_name_rm)
        print("upload", file_name_rm)
        command = "rm -rf " + file_name_rm
        os.system(command)
    except:
        continue