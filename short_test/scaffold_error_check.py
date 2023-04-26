import os
from os import listdir
from os.path import isfile, join
import time
import pandas as pd
import numpy as np
import multiprocessing

from datetime import date
from multiprocessing import Process, current_process

from rdkit import Chem, RDConfig
from rdkit.Chem import MolStandardize, Draw, AllChem, rdMolAlign
from rdkit.Chem.Scaffolds import rdScaffoldNetwork
from Bio.PDB import PDBParser, PDBIO, Select
import Bio.PDB
from google.cloud import storage
import time
import psycopg2
from sqlalchemy import desc, select
import sqlalchemy
import re


smi =" CCC1O[C@H](O[Si](C)(C)C(C)(C)C)C(NS(=O)(=O)c2ccccc2)[C@H](OCc2ccccc2)[C@@H]1O[C@@H]1OC(COCc2ccccc2)[C@@H](O[C@@H]2OC(CO[C@H]3OC(COCc4ccccc4)[C@@H](OCc4ccccc4)[C@@H](OCc4ccccc4)C3OC3CC(COCc4ccccc4)[C@@H](OC4CC(COCc5ccccc5)[C@H](OC(C)=O)[C@H](O[C@]5(C(=O)OC)C[C@@H](OC(C)=O)[C@@H](NC(C)=O)C([C@H](OC(C)=O)[C@@H](COC(C)=O)OC(C)=O)O5)[C@@H]4OCc4ccccc4)[C@H](OCc4ccccc4)[C@@H]3C)[C@@H](OCc3ccccc3)[C@H](OC3O[C@@H](COCc4ccccc4)[C@@H](OCc4ccccc4)C(OCc4ccccc4)C3OC3OC(COCc4ccccc4)[C@@H](OC4CC(COCc5ccccc5)[C@H](OC(C)=O)[C@H](O[C@]5(C(=O)OC)C[C@@H](OC(C)=O)[C@@H](NC(C)=O)C([C@H](OC(C)=O)[C@@H](COC(C)=O)OC(C)=O)O5)[C@@H]4OCc4ccccc4)[C@H](OCc4ccccc4)[C@@H]3C)C2OCc2ccccc2)[C@@H](OCc2ccccc2)C1NS(=O)(=O)c1ccccc1"
print('smi ok')
mol = Chem.MolFromSmiles(smi)
print('mol ok')
params = rdScaffoldNetwork.ScaffoldNetworkParams()
params.includeScaffoldsWithAttachments=False

try:
    netwks = rdScaffoldNetwork.CreateScaffoldNetwork([mol], params)
except:
    print(" not netwks")


node_list = []
edge_list = []
frag_list = {}

for edge in netwks.edges:
    if edge.beginIdx != 0:
        edge_list.append([edge.beginIdx - 1, edge.endIdx - 1 ])

for i, node in enumerate(netwks.nodes):
    if i != 0:
        node_list.append([node, netwks.edges[i-1].type])

print(node_list)