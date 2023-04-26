from time import perf_counter
import psycopg2
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import MolStandardize
from sqlalchemy import create_engine, engine
from sqlalchemy.orm import sessionmaker
import pandas as pd

jpc_info = "host='34.64.203.48' dbname='chemdb' user='postgres' password='1234' port='5433'"
conn_jpc = psycopg2.connect(jpc_info)
cur_jpc = conn_jpc.cursor()

smi = 'Cn1ccnc1Sc1cc(C(=O)N=c2[nH]ccs2)c(N)cc1F'
percentage = str(0.8)

command = "SELECT *, canonical_smiles |~| '" + smi + "' FROM mirror_mol2 WHERE ('" + smi + "'," + percentage + ")::sim_order |<~| canonical_smiles;" 


cur_jpc.execute(command)
print('connect')
result_frag = cur_jpc.fetchall()

id_list = []
smile_list = []

for result in result_frag:
    print(result[0], result[1])
    id_list.append(result[0])
    smile_list.append(result[1])

df = pd.DataFrame()
df['id']=id_list
df['smiles']=smile_list
df.to_csv('/Users/song-inhyeok/Documents/SCH_API_PROJECT/1V4S_smi08.csv', index = False)    
    
print('end')