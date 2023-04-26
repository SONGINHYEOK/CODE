import requests
URL = "https://zinc20.docking.org/substances.smi?count=all"
  
file = requests.get(URL, stream = True, timeout=None)


with requests.Session() as s:
    r = s.get(URL)
    
    with open("./ZINC20_ALL.smi","wb") as pdf:
        try:
            for chunk in r.iter_content(chunk_size=200000):
            
                if chunk:
                    pdf.write(chunk)
        except:
            pass        
            
  
