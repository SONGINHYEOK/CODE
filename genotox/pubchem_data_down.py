from bs4 import BeautifulSoup as bs
import requests
from urllib.request import urlretrieve

fix_path = "https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/CSV/Description/"

response = requests.get(fix_path)
html_txt = response.text

soup = bs(response.text, 'html.parser')

href_sel = soup.select("a")

href_list = []

for i in href_sel:
    href  = i.attrs['href']
    href_list.append(href)

del href_list[0]
del href_list[-1]




import os

# path_folder의 경로는 각자 저장할 폴더의 경로를 적어줄 것(ex.img_download)
path_folder = '/Users/song-inhyeog/CODEING/CODE/genotox/data/'


    
#link 만들어서 다운 받기

down_link = []

for hr in href_list:
    fi_link = fix_path + hr
    down_link.append(fi_link)


m = 0
for dw in down_link:
    urlretrieve(dw, path_folder + dw.split('/')[-1])
    m +=1
    print(m)