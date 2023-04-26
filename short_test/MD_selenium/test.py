from selenium import webdriver
from selenium.webdriver.chrome.options import Options

from webdriver_manager.chrome import ChromeDriverManager
import time

from bs4 import BeautifulSoup
import os, sys, copy

import shutil
from glob import glob
import os
from os import listdir
from os.path import isfile, join
import time
import paramiko
import getpass
from scp import SCPClient

options = webdriver.ChromeOptions()

options.add_argument('window-size=1920x1080')



driver = webdriver.Chrome(executable_path=ChromeDriverManager().install(), chrome_options=options)


driver.implicitly_wait(0.5)
driver.maximize_window()

driver.get("http://34.64.245.194:9000")

driver.find_element_by_xpath('//*[@id="app"]/div/div/div[1]/div[3]/div/div/div[1]/button[1]').click()

time.sleep(1)

driver.find_element_by_xpath('//*[@id="app"]/div/div/div[1]/div[3]/div/div/div[2]/div[1]/div[5]/div/button').click()
time.sleep(1)


pdb = driver.find_element_by_xpath('//*[@id="app"]/div/div/div[1]/div[3]/div/div/div[2]/div[1]/div[5]/div[3]/div/div/input')
pdb.send_keys("/Users/song-inhyeok/CODING/short_test/1V4S_MD_result.pdb")

time.sleep(1)

xtc = driver.find_element_by_xpath('//*[@id="app"]/div/div/div[1]/div[3]/div/div/div[2]/div[1]/div[5]/div[4]/div/div/input')
xtc.send_keys("/Users/song-inhyeok/CODING/short_test/1V4S_MD_result.xtc")

time.sleep(1)

driver.find_element_by_xpath('//*[@id="app"]/div/div/div[1]/div[3]/div/div/div[2]/div[1]/div[5]/div[5]/div/button').click()

time.sleep(1)


#wremove water ligand focus
driver.find_element_by_xpath('//*[@id="app"]/div/div/div[1]/div[3]/div/div/div[1]/button[7]').click()
time.sleep(1)
driver.find_element_by_xpath('//*[@id="app"]/div/div/div[1]/div[4]/div/div/div[4]/div[1]/button').click()
time.sleep(1)
driver.find_element_by_xpath('//*[@id="app"]/div/div/div[1]/div[4]/div/div/div[4]/div[3]/div[2]/button[1]').click()
time.sleep(1)
driver.find_element_by_xpath('//*[@id="app"]/div/div/div[1]/div[4]/div/div/div[4]/div[3]/div[3]/button[2]').click()
time.sleep(1)

#upload 

driver.find_element_by_xpath('//*[@id="app"]/div/div/div[1]/div[1]/div/div[6]/button').click()

driver.find_element_by_xpath('//*[@id="app"]/div/div/div[1]/div[6]/div/div/div[2]/div/button').click()

server_url=driver.find_element_by_xpath('//*[@id="app"]/div/div/div[1]/div[6]/div/div/div[2]/div[2]/div/input')
server_url.clear()
name=driver.find_element_by_xpath('//*[@id="app"]/div/div/div[1]/div[6]/div/div/div[2]/div[3]/div/input')
description=driver.find_element_by_xpath('//*[@id="app"]/div/div/div[1]/div[6]/div/div/div[2]/div[4]/div/input')
source=driver.find_element_by_xpath('//*[@id="app"]/div/div/div[1]/div[6]/div/div/div[2]/div[5]/div/input')

server_url.send_keys("http://34.64.245.194:9001")

job_id = 'focus_test'
name.send_keys(job_id)
description.send_keys(job_id)
source.send_keys(job_id)

time.sleep(1)

driver.find_element_by_xpath('//*[@id="app"]/div/div/div[1]/div[6]/div/div/div[2]/div[6]/button').click()

time.sleep(5)


"""
session_url =driver.find_element_by_xpath('//*[@id="app"]/div/div/div[1]/div[3]/div/div/div[2]/div[3]/div/div[2]/div[1]/div/input')
session_url.send_keys("http://34.64.245.194:9001")
time.sleep(2)


driver.find_element_by_xpath('//*[@id="app"]/div/div/div[1]/div[3]/div/div/div[2]/div[3]/div/div[2]/div[2]/button').click()

html =  driver.page_source
soup = BeautifulSoup(html,'html.parser')

print(soup)
"""
#search_list = soup.find_all('li')
#for i in search_list:
 #   print(i)



#my_element = driver.find_elements_by_xpath(f"//button[contains(text(),{job_id})]")
#my_element.send_keys('data=id')
#my_element = driver.find_element_by_link_text(job_id)

#print(my_element)
#print(my_element.get_attribute('data-id'))




#driver.close()



#to identify element
#s = driver.find_element_by_xpath("//input[@type='file']")
#file path specified with send_keys
#s.send_keys("C:\Users\Pictures\Logo.jpg")