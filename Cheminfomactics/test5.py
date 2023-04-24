import urllib.request
import requests
#'https://zinc20.docking.org/substances/ZINC002326992000.smi'

import pandas as pd
import os
from selenium import webdriver

options = webdriver.ChromeOptions()
options.add_argument('headless')
options.add_argument('window-size=1920x1080')
options.add_argument("disable-gpu")
options.add_experimental_option("prefs", {
  "download.default_directory": r"/Users/song-inhyeok/Documents/ZINC20_cr_project/ZINC_ALL/",
  "download.prompt_for_download": False,
  "download.directory_upgrade": True,
  "safebrowsing.enabled": True
})

path = '/Users/song-inhyeok/Documents/coding/PDB/chromedriver'    
save_path = "/Users/song-inhyeok/Documents/data/ZINC20_data/" 

root_url = 'https://zinc20.docking.org/substances/'
root_file_path = '/Users/song-inhyeok/Documents/ZINC20_cr_project/ZINC20_smi/18/'
file_list=os.listdir(root_file_path)
#print(file_list)



driver = webdriver.Chrome(path, chrome_options=options)

count = 0
for file in file_list:
    df = pd.read_csv(root_file_path+file, header=False)
    df

    #print(zid)
    """
    for i in zid:
        url = root_url + i           
        url = str(url)
        print(type(url))
        driver.get(url)
        driver.implicitly_wait(0.5)
    count +=1

    print(count)
    """

