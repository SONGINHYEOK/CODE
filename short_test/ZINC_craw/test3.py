import requests

session = requests.Session()
session.post('https://zinc20.docking.org/substances.smi?count=all')



file_url = "https://zinc20.docking.org/substances.smi?count=all"
resp = session.get(file_url, stream=True)

if resp.status_code == 200:
    with open("./ZINC20_ALL.smi","wb") as pdf:
        try:
            for chunk in resp.iter_content(chunk_size=200000):
            
                if chunk:
                    pdf.write(chunk)
        except:
            pass      