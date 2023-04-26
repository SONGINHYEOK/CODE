
import os
from this import d
import time
import argparse
from matplotlib.pyplot import axis
from sqlalchemy import create_engine, engine
from sqlalchemy.orm import sessionmaker
from dbmodel import ConformerFingerprint
import csv
import pandas as pd
import psycopg2

db = create_engine(
        engine.url.URL(
            drivername='postgresql',
            username='chemdb-admin',
            password='net-targets:0813',
            #host='34.64.95.244', # Chemdb (protein), dynamichemdb 
            host='34.64.79.156', # Chemdb -  services DB /molecule 
            #host='35.193.86.103', # Chemdb - backup
            port='5432',
            #database='dynamichemdb' # fragment
            database='uptest' # Chemdb - protein  / molecule
            #database='frag_test'
        ),
    )

db.echo = False
SessionFactory = sessionmaker(bind=db, autoflush=False)
session = SessionFactory()


start = time.time()

with open("/Users/song-inhyeok/Documents/coding/sample_40.csv", 'r') as file:
    data_df = pd.read_csv(file)
    data_df.head()
#data_df.to_sql('conformer_fingerprint2', con=db, index=False, if_exists='append')
end = time.time() - start



    
print('file read Running time;', end)

\copy conformer_fingerprint(fp, aa, ad, ah, an, ap, ar, dd, dh, dn, dp, dr, hh, hn, hp, hr, nn, np, nr, pp, pr, rr, max_elements) FROM '/Users/song-inhyeok/Documents/coding/sample_test.csv' DELIMITER ',' encoding 'utf8';


\copy conformer_fingerprint2(fp, aa, ad, ah, an, ap, ar, dd, dh, dn, dp, dr, hh, hn, hp, hr, nn, np, nr, pp, pr, rr, max_elements) FROM '/Users/song-inhyeok/Documents/coding/sample_test.csv' DELIMITER ',' encoding 'utf8';