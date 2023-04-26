from audioop import rms
from logging import root
import os
from os import listdir, system
from os.path import isfile, join
import re
import readline
from symbol import sym_name
from numpy import source

import rdkit
import pytz
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem, MolStandardize, CombineMols
from rdkit.Chem.Scaffolds import rdScaffoldNetwork

import copy
import pandas as pd
import time
from datetime import datetime

import Bio.PDB
import psycopg2

from Bio.PDB import PDBParser, PDBIO, Select
import paramiko
import getpass
from scp import SCPClient
import tarfile
from sqlalchemy import desc, select
import sqlalchemy
import shutil

score_list =[]
id_list = []
info_dict = {}
root_path = "/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/docking"
file_list=os.listdir(root_path)

for file in file_list:
    file_name = root_path + '/' + file
    with open(file_name, 'r') as fi:
        line = fi.readline().rstrip()
        id_list.append(line)
    sdf = Chem.SDMolSupplier(file_name)
    #print(sdf['r_i_docking_score'])
    for i in sdf:
        score_list.append(i.GetProp('r_i_docking_score'))

for id, score in zip(id_list, score_list):
    info_dict[id]=score

print(info_dict)