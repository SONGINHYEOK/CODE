
import os 
import sys
from turtle import pd
from prody import *
from rdkit import Chem
from rdkit.Chem import AllChem
from io import StringIO
import pypdb
from Bio.PDB import PDBParser, PDBIO, Select
import pandas as pd

path = '/Users/song-inhyeok/Documents/PDB/lignad'


for pfb_file in os.listdir(path):
    i = 1
    if pfb_file.endswith('.pdb') and not pfb_file.startswith("lig_"):
        pdb_code = pfb_file[:-4]
        pdb = parsePDB(pdb_code)
        for model in pdb:
            model_atoms = Bio.PDB.Selection.unfold_entities(model, 'A')
            for chain in model:
                for residue in chain:
                    start_lig = time.time()
                    if not is_het(residue):
                        continue
        
        
        output = StringIO()
        sub_mol = ligand.select(f"resname {res_name}")
        chem_desc = pypdb.describe_chemical(f"{res_name}")
        sub_smiles = chem_desc["describeHet"]["ligandInfo"]["ligand"]["smiles"]
        template = AllChem.MolFromSmiles(sub_smiles)
        writePDBStream(output, sub_mol)
        pdb_string = output.getvalue()
        rd_mol = AllChem.MolFromPDBBlock(pdb_string)
        new_mol = AllChem.AssignBondOrdersFromTemplate(template, rd_mol)
        
        
        
        
output = StringIO()
sub_mol = ligand.select(f"resname {res_name}")
chem_desc = pypdb.describe_chemical(f"{res_name}")
"""        
  
        
def get_ligands(pdb_id):
         """Return ligands of given PDB ID

     Parameters
     ----------

     pdb_id : string
         A 4 character string giving a pdb entry of interest

     Returns
     -------

     out : dict
         A dictionary containing a list of ligands associated with the entry

     Examples
     --------
     >>> ligand_dict = get_ligands('100D')
     >>> print(ligand_dict)
     {'id': '100D',
     'ligandInfo': {'ligand': {'@chemicalID': 'SPM',
                            '@molecularWeight': '202.34',
                            '@structureId': '100D',
                            '@type': 'non-polymer',
                            'InChI': 'InChI=1S/C10H26N4/c11-5-3-9-13-7-1-2-8-14-10-4-6-12/h13-14H,1-12H2',
                            'InChIKey': 'PFNFFQXMRSDOHW-UHFFFAOYSA-N',
                            'chemicalName': 'SPERMINE',
                            'formula': 'C10 H26 N4',
                           'smiles': 'C(CCNCCCN)CNCCCN'}}}

     """
    out = get_info(pdb_id, url_root = 'http://www.rcsb.org/pdb/rest/ligandInfo?structureId=')
#     out = to_dict(out)
#     return remove_at_sign(out['structureId']

from collections import OrderedDict, Counter
from itertools import repeat, chain
import requests
import time
import re
import json
import warnings

from pypdb.util import http_requests

def to_dict(odict):
    '''Convert OrderedDict to dict

    Takes a nested, OrderedDict() object and outputs a
    normal dictionary of the lowest-level key:val pairs

    Parameters
    ----------

    odict : OrderedDict

    Returns
    -------

    out : dict

        A dictionary corresponding to the flattened form of
        the input OrderedDict

    '''

    out = json.loads(json.dumps(odict))
    return out

from pandas import json_normalize


def get_info(pdb_id, url_root='https://data.rcsb.org/rest/v1/core/entry/'):
    '''Look up all information about a given PDB ID

    Parameters
    ----------

    pdb_id : string
        A 4 character string giving a pdb entry of interest

    url_root : string
        The string root of the specific url for the request type

    Returns
    -------

    out : dict()
        An ordered dictionary object corresponding to entry information

    '''
    pdb_id = pdb_id.replace(":", "/")  # replace old entry identifier
    url = url_root + pdb_id
    response = http_requests.request_limited(url)

    if response is None or response.status_code != 200:
        warnings.warn("Retrieval failed, returning None")
        return None

    result = str(response.text)

    out = json.loads(result)

    return out



get_info('6yp7')



get_all_info = get_info  # Alias
describe_pdb = get_info  # Alias for now; eventually make this point to the Graph search https://data.rcsb.org/migration-guide.html#pdb-file-description
get_entity_info = get_info  # Alias


def get_ligands(pdb_id):
    """Return ligands of given PDB ID

     Parameters
     ----------

     pdb_id : string
         A 4 character string giving a pdb entry of interest

     Returns
     -------

     out : dict
         A dictionary containing a list of ligands associated with the entry

     Examples
     --------
     >>> ligand_dict = get_ligands('100D')
     >>> print(ligand_dict)
     {'id': '100D',
     'ligandInfo': {'ligand': {'@chemicalID': 'SPM',
                            '@molecularWeight': '202.34',
                            '@structureId': '100D',
                            '@type': 'non-polymer',
                            'InChI': 'InChI=1S/C10H26N4/c11-5-3-9-13-7-1-2-8-14-10-4-6-12/h13-14H,1-12H2',
                            'InChIKey': 'PFNFFQXMRSDOHW-UHFFFAOYSA-N',
                            'chemicalName': 'SPERMINE',
                            'formula': 'C10 H26 N4',
                           'smiles': 'C(CCNCCCN)CNCCCN'}}}

    """
    out = get_info(pdb_id, url_root = 'http://www.rcsb.org/pdb/rest/ligandInfo?structureId=')
    out = to_dict(out)

get_info('6yp7')
get_ligands('6yp7')