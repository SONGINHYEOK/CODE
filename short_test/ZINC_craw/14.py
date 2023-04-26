import urllib.request
import requests
#'https://zinc20.docking.org/substances/ZINC002326992000.smi'


from selenium import webdriver

options = webdriver.ChromeOptions()
options.add_argument('headless')
options.add_argument('window-size=1920x1080')
options.add_argument("disable-gpu")
options.add_experimental_option("prefs", {
  "download.default_directory": r"/Users/song-inhyeok/Documents/data/ZINC20_data/",
  "download.prompt_for_download": False,
  "download.directory_upgrade": True,
  "safebrowsing.enabled": True
})

path = '/Users/song-inhyeok/Documents/coding/PDB/chromedriver'    
save_path = "/Users/song-inhyeok/Documents/data/ZINC20_data/" 

def list_chunk(lst, n):
    return [lst[i:i+n] for i in range(0, len(lst), n)]


list_test = list(range(1,23269921))
list_chunked = list_chunk(list_test, 1163496)  

count = 0
driver = webdriver.Chrome(path, chrome_options=options)

link_list = []
for i in list_chunked[13]:
    url='https://zinc20.docking.org/substances/ZINC'+str(i).zfill(12)+'.smi'
    driver.get(url)
    driver.implicitly_wait(0.5)
    count +=1

    print(count)


