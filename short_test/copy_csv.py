from base64 import encode
from email.header import Header
import imp
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
            database='uptest', # Chemdb - protein  / molecule
            #database='frag_test'
            charset='utf8'
        ),
    )

db.echo = False
SessionFactory = sessionmaker(bind=db, autoflush=False)
session = SessionFactory()
conn = psycopg2.connect(db)
cur = conn.cursor()

from sqlalchemy import text

command = "copy conformer_fingerprint2(fp, aa, ad, ah, an, ap, ar, dd, dh, dn, dp, dr, hh, hn, hp, hr, nn, np, nr, pp, pr, rr, max_elements) FROM '\Users\song-inhyeok\Documents\coding\sample_40.csv' DELIMITER ',' encoding 'utf8';
cur.execute(command)
conn.commit()
cur.close()
conn.close()



"""
f = open('/Users/song-inhyeok/Documents/coding/sample_40.csv', 'r', encoding='utf-8')
#rdr = csv.reader(f)
rdr = f.readlines()
"""


