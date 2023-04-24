from asyncore import write
from concurrent.futures import process
from glob import glob
import os
from os import listdir
from os.path import isfile, join
import time
import pandas as pd
import numpy as np
import multiprocessing
from multiprocessing import Process, current_process
from schprocess import *
from datetime import date
from featureutil import make_pair_distance, make_distance_fingerprint, parse_feature_sch
from extractprotein import *
from standardizeutil import *
from Bio.PDB import PDBParser, PDBIO, Select
import Bio.PDB
from rdkit import Chem, RDConfig
from rdkit.Chem import MolStandardize, Draw, AllChem, rdMolAlign
import psycopg2
from sqlalchemy import desc, select
import sqlalchemy
from google.cloud import storage
from openbabel import openbabel
import os, sys, copy
from openbabel import pybel
import subprocess
from dbmodel import Protein_structure, Ligand, Binding_site, Site_residue, Amino_acid, Molecule, Result_conformer, Job, Conformer
import time
import shutil
from datetime import datetime
import pytz

class NonHetSelect(Select):
    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0

class Make_system:
    save_dir = '/home/NCAP_test/media'
    KST = 'Asia/Seoul'
    read_path = "PDB/"
    #read_path2 = "PDB_prep/"
    save_path = "./"
    
    def __init__(self, session, file_name, grid_size, sol):#, bsite_id):#, , job_id):
        self.session = session
        self.file_name = file_name
        self.grid_size = grid_size
        self.sol = sol
    def act_system(self):
    
        protein_bucket_name = 'nettargets-discoveryweb-protein-bucket'
        os.environ["GCLOUD_PROJECT"] = "nettargets-chemdb-v1"
        client = storage.Client()
        bucket = client.get_bucket(protein_bucket_name)
        
        for pfb_file in bucket.list_blobs(prefix=self.read_path + self.file_name):
            if pfb_file.name.endswith('_prepped.pdb') == False:
                pdb_file_name = pfb_file.name.replace(self.read_path, "")
                print(pdb_file_name)
                pdb_code = pdb_file_name.replace(".pdb","")
                save_pdb = self.save_path + pfb_file.name.replace('PDB/', '')
                print(save_pdb, pfb_file.name)
                pfb_file.download_to_filename(save_pdb)
    
        #addH_pdb = self.save_path +f"{self.file_name}_rec.pdb"
        #bs = self.session.query(Binding_site).filter(Binding_site.id==self.bsite_id).first() 
        #li=self.session.query(Ligand).filter(Ligand.ligand_id==bs.ligand.ligand_id).filter(Ligand.id==bs.ligand.id).first() 
        
        
        procsch_addh_rec(save_pdb, save_pdb)
        
        pdb = PDBParser().get_structure(pdb_code, save_pdb)
        io = PDBIO()
        io.set_structure(pdb)
        
        for model in pdb:
            io.save(self.save_path+f"/{self.file_name}_addH_rec.pdb", NonHetSelect())
          
        target_mol = "target.sdf"    
        root = "/home/NCAP_test/"
        chain_file = self.save_path+f"/{self.file_name}_addH_rec.pdb"
       
          
        sdf_path = root=target_mol
        tmp_addH_file = self.save_path+f"{self.file_name}_addH_rec_merge.pdb"
        
        merge_file(chain_file, sdf_path, tmp_addH_file)
        
        replace_in_file(tmp_addH_file, "UNK", "LIG")
        
        tmp_addH_file = self.save_path+f"{self.file_name}_addH_rec_merge.pdb"
        tmp_addH_rech = f"{self.file_name}_addH_rec_merge.maegz"
        procsch_addh_rec(tmp_addH_file, tmp_addH_rech)
    
        output_name = pdb_code +'_prep.maegz'
        
        #protein preparation
        protein_prep(tmp_addH_rech, output_name)
        while os.path.exists('/home/NCAP_test/'+f"{self.file_name}_prep.maegz") !=True:
            time.sleep(10)
            if os.path.exists('/home/NCAP_test/'+f"{self.file_name}_prep.maegz") ==True:
                break
        
        print('complete prep')
        
        with open(f"{self.file_name}.msj", 'w') as msj:
            msj.write("task {"  "\n")
            msj.write("   task = 'desmond:auto'" "\n")
            msj.write("}" "\n")

            msj.write("build_geometry {" "\n")
            msj.write("    add_counterion = {" "\n")
            msj.write("            ion = Na" "\n")
            msj.write("           number = neutralize_system" "\n")
            msj.write("}          " "\n")
            msj.write("   box = {" "\n")
            msj.write("       shape = orthorhombic" "\n")
            msj.write(f"       size = [{self.grid_size} {self.grid_size} {self.grid_size} ]" "\n")
            msj.write("       size_type = buffer " "\n")
            msj.write("}               " "\n")
            msj.write('    override_forcefield = S-OPLS'"\n")
            msj.write("rezero_system = true" "\n")
            msj.write(f"solvent = {self.sol}" "\n")
            msj.write("}"               "\n"  )
            msj.write('assign_forcefield {' "\n")        
            msj.write("                  forcefield = S-OPLS" "\n")
            msj.write("}")
            msj.close()

        time.sleep(5)
       
        make_model(self.file_name)    

        while os.path.exists('/home/NCAP_test/'+f"{self.file_name}-out.cms") !=True:
            time.sleep(10)
            if os.path.exists('/home/NCAP_test/'+f"{self.file_name}-out.cms") ==True:
                break
        print('complete system')
        
def replace_in_file(file_path, old_str, new_str):
        # 파일 읽어들이기
        fr = open(file_path, 'r')
        lines = fr.readlines()
        fr.close()
        
        # old_str -> new_str 치환
        fw = open(file_path, 'w')
        for line in lines:
            fw.write(line.replace(old_str, new_str))
        fw.close()
