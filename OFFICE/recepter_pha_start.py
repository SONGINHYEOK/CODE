from ast import Num
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
from dbmodel import Protein_structure, Ligand, Binding_site, Site_residue, Amino_acid, Molecule

"""
#======================== Workflow ================================#
1.gcp storage get pdb file
2.schrodinger input pdb file ->  output file receptor.mae (only protein)
3.db search protein_pdb -> protein_structure -> binding_site -> site_residue get residue_number
4.ressidue_number -> make "schrodinger ASL" -> get -site_center
5.schrodinger make receptor pharmacophore process -> out recpharm.phypo
#=================================================================#
"""
class NonHetSelect(Select):
    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0

class Receptor_proc:
    save_dir = './BIChem/media/'
    KST = 'Asia/Seoul'
    read_path = "PDB/"
    save_path = "./"
    
    def __init__(self, session, file_name, bsite_id):
        self.session = session
        self.file_name = file_name
        self.bsite_id = bsite_id
    def receptor_pharm(self):
        num_list = []
        
        #download pdb file
        
              
        for pfb_file in bucket.list_blobs(prefix=self.read_path + self.file_name):
            pdb_file_name = pfb_file.name.replace(self.read_path, "")
            pdb_code = pdb_file_name.replace(".pdb","")
            save_pdb = self.save_path + pfb_file.name.replace('PDB/', '')
            print(save_pdb, pfb_file.name)
            pfb_file.download_to_filename(save_pdb)
        
        #addH_pdb = self.save_path +f"{pdb_code}_rec.pdb"
        procsch_addh_rec(save_pdb, save_pdb)
        i = 1
        pdb = PDBParser().get_structure(pdb_code, save_pdb)
        io = PDBIO()
        io.set_structure(pdb)
        for model in pdb:
            io.save(self.save_path+f"/{pdb_code}_addH_rec.pdb", NonHetSelect())
            
            for chain in model:
                for residue in chain:
                        original_bool = False
                        if not is_het(residue):
                            continue
                        chain_id = chain.id[0]
                        residue_id = residue.id[1]
                        #print("save", f"lig_{pdb_id}_{residue.id[0]}_{residue.id[1]}_{i}")
                        ligand_id = residue.id[0].split('_')[1]
                        print(chain_id, residue_id, ligand_id)
                        #io.save(f"lig_{pdb_id}_{residue.id[0]}_{i}.pdb", ResidueSelect(chain, residue))
                        io.save(self.save_path+f"lig_{pdb_code}_{chain_id}_{ligand_id}_{residue_id}_{i}.pdb", ResidueSelect(chain, residue))

                        obConversion = openbabel.OBConversion()
                        obConversion.SetInAndOutFormats("pdb", "sdf")

                        mol = openbabel.OBMol()

                        li_file_name = f"lig_{pdb_code}_{chain_id}_{ligand_id}_{residue_id}_{i}"

                        obConversion.ReadFile(mol, self.save_path+li_file_name+".pdb"  ) # Open Babel will uncompress automatically
                        mol.AddHydrogens()

                        obConversion.WriteFile(mol,self.save_path+li_file_name+".sdf")        
                        os.remove(self.save_path+li_file_name+".pdb")
            
        
        pt_li_bs = self.session.query(Binding_site).filter(Binding_site.id==self.bsite_id).first()
                
        tmp_lig_file = self.save_path+f"lig_{pdb_code}_{chain_id}_{pt_li_bs.ligand.ligand_id}_{pt_li_bs.ligand.residue_number}_{i}"+".sdf"
        tmp_lig_maegz = self.save_path+f"lig_{pdb_code}_{chain_id}_{pt_li_bs.ligand.ligand_id}_{pt_li_bs.ligand.residue_number}_{i}"+".maegz"
        procsch_maegztosdf(tmp_lig_file, tmp_lig_maegz)
         
        tmp_addH_file = self.save_path+f"{pdb_code}_addH_rec.pdb"
        tmp_addH_rech = f"{pdb_code}_addH_rec.maegz"
        procsch_addh_rec(tmp_addH_file, tmp_addH_rech)
                    
        #make phypo
        #procsch_recpharm(tmp_addH_rech, tmp_lig_maegz, center,pdb_code)
        procsch_recpharm(tmp_addH_rech, tmp_lig_maegz ,pdb_code)
        time.sleep(10)
        
        #phypo convert txt
        feat_txt=procsch_vsquery(f'./Hypothesis_{pdb_code}_receptor.phypo')
        pharm_file = open('recp_pharm.txt', 'w')
        pharm_file.write(feat_txt +'\n')
        pharm_file.close()
        
        #remove tmp make file
        os.remove(save_pdb)
        os.remove(tmp_addH_rech)
        os.remove(tmp_addH_file)
        os.remove(tmp_lig_maegz)
        self.remove_tmp_sdf('./', '.sdf')
        
        os.remove(f"Hypothesis_{pdb_code}.phypo")
        os.remove(f"Hypothesis_{pdb_code}_receptor.phypo")
        
    def remove_tmp_sdf(self,filePath, fileExtension):
        if os.path.exists(filePath):
            for file in os.scandir(filePath):
                if file.name.endswith(fileExtension):
                    os.remove(file.path)
    
