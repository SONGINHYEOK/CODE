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

class glide_proc_test:
    save_dir = '/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/media'
    KST = 'Asia/Seoul'
    read_path = "PDB/"
    #read_path2 = "PDB_prep/"
    save_path = "./"
    
    def __init__(self, session, file_name, bsite_id, job_id, grid_size, grid_center, precis):
        self.session = session
        self.file_name = file_name
        self.bsite_id = bsite_id
        self.job_id = job_id
        self.grid_size = grid_size
        self.grid_center = grid_center
        self.precis = precis
    
    def finish_job(self):
        # self.session.rollback() # temp process - query feature included in session
        KST = pytz.timezone(self.KST)
        end_dt = datetime.now(KST)
        job = self.session.query(Job).filter_by(id=self.job_id).first()
        job.end_date = end_dt
        self.session.commit()
        self.session.close()
    
    def glide_process_test(self):

        num_list = []
       
        #download pdb file
       
        """
        for pfbp_file in bucket.list_blobs(prefix=self.read_path2 + self.file_name):
            pdbp_file_name = pfbp_file.name.replace(self.read_path2, "")
            print(pdbp_file_name)
                 
            save_pdbp = self.save_path + pfbp_file.name.replace('PDB_prep/', '')
            print(save_pdbp, pfbp_file.name)
            pfbp_file.download_to_filename(save_pdbp)
        """
        for pfb_file in bucket.list_blobs(prefix=self.read_path + self.file_name):
            if pfb_file.name.endswith('_prepped.pdb') == False:
                pdb_file_name = pfb_file.name.replace(self.read_path, "")
                print(pdb_file_name)
                pdb_code = pdb_file_name.replace(".pdb","")
                save_pdb = self.save_path + pfb_file.name.replace('PDB/', '')
                print(save_pdb, pfb_file.name)
                pfb_file.download_to_filename(save_pdb)
    
        #addH_pdb = self.save_path +f"{self.file_name}_rec.pdb"
        procsch_addh_rec(save_pdb, save_pdb)
        
        pdb = PDBParser().get_structure(pdb_code, save_pdb)
        io = PDBIO()
        io.set_structure(pdb)
        for model in pdb:
            io.save(self.save_path+f"/{self.file_name}_addH_rec.pdb", NonHetSelect())
        
        prt = self.session.query(Protein_structure).filter(Protein_structure.pdb_id==self.file_name).first()
        bs = self.session.query(Binding_site).filter(Binding_site.id==self.bsite_id).filter(Binding_site.protein_str_id==prt.id).first()
        sites=bs.site_residue       
        li=self.session.query(Ligand).filter(Ligand.ligand_id==bs.ligand.ligand_id).filter(Ligand.id==bs.ligand_id).first() 

        tmp_addH_file = self.save_path+f"{self.file_name}_addH_rec.pdb"
        tmp_addH_rech = f"{self.file_name}_addH_rec.maegz"
        procsch_addh_rec(tmp_addH_file, tmp_addH_rech)
        
        
    
        output_name = pdb_code +'_prep.maegz'
        
        #protein preparation
        protein_prep(save_pdb, output_name)
        while os.path.exists('/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/'+f"{self.file_name}_prep.maegz") !=True:
            time.sleep(10)
            if os.path.exists('/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/'+f"{self.file_name}_prep.maegz") ==True:
                break
        
        print('complete prep')
        if self.grid_center == 'Ligand':
            #make res_asl 
            res_asl = f"res.num {li.residue_number}" 
            print(res_asl)
            chain_asl = f"chain.name  {li.chain}"
            final_asl = f"({res_asl}) AND (({chain_asl}))"
            center=procsch_recsite(output_name, final_asl)
            print(center)
    
        else:    
            for site in sites:    
                num_list.append(site.residue_number)
        
            res_asl = "res.num " + ','.join(str(e) for e in num_list)
            print(res_asl)
            chain_asl = f"chain.name  {site.chain}"
            final_asl = f"({res_asl}) AND (({chain_asl}))"
            center=procsch_recsite(output_name, final_asl)
            print(center)
        
      
        #make grid in file  
        with open(f"./{self.file_name}.in", 'a') as f:
            f.write('FORCEFIELD   OPLS_2005' +'\n')
            f.write(f'GRID_CENTER   {center}') 
            f.write(f'GRIDFILE   {self.file_name}_grid.zip'+'\n')
            f.write('INNERBOX   10, 10, 10'+'\n')    
            f.write(f'OUTERBOX   {self.grid_size}, {self.grid_size}, {self.grid_size}'+'\n')
            f.write(f'RECEP_FILE   {output_name}')
            f.close()
        
        while os.path.exists('/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/'+f"{self.file_name}.in") !=True:
            time.sleep(10)
            if os.path.exists('/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/'+f"{self.file_name}.in") ==True:
                break
        
        print('complete make in file')
        
        
        #make grid
        
        make_grid(f"./{self.file_name}.in")
        
        while os.path.exists('/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/'+f"{self.file_name}_grid.zip") !=True:
            time.sleep(10)
            if os.path.exists('/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/'+f"{self.file_name}_grid.zip") ==True:
                break
        
        print('complete make grid')
        
        #make glide docking in file
        with open(f"./{self.file_name}_docking.in", 'a') as doc:
            doc.write('FORCEFIELD   OPLS_2005'+'\n')
            doc.write('GRIDFILE' +'  '+f'/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/{self.file_name}_grid.zip'+'\n')
            doc.write('LIGANDFILE   /Users/song-inhyeok/CODING/DRUGBANK/3D_structures.sdf'+'\n')
            doc.write(f'PRECISION   {self.precis}')
            doc.close()
        
        print('start glide docking')
        glide_start(self.file_name+"_docking.in")
        
        while os.path.exists('/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/'+self.file_name+"_docking_pv.maegz") !=True:
            time.sleep(10)
            if os.path.exists('/Users/song-inhyeok/CODING/PROTOTYPE/BIChem'+self.file_name+"_docking_pv.maegz") ==True :
                break
        print('end glide docking')
        
        #split hit maegz - > sdf
        procsch_pvsdf(f"{self.file_name}_docking_pv.maegz", self.file_name+"_docking.sdf")    
        
        #result save
        score_list =[]
        id_list = []
        info_dict = {}
        root_path = "/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/docking"
        file_list=os.listdir(root_path)

        for file in file_list:
            file_name = root_path + '/' + file
            with open(file_name, 'r') as fi:
                line = fi.readline().rstrip()
            sdf = Chem.SDMolSupplier(file_name)
            #print(sdf['r_i_docking_score'])
            for i in sdf:
                #score_list.append(i.GetProp('r_i_docking_score'))
                glide_score = i.GetProp('r_i_docking_score')
                mol_3d_block= Chem.MolToMolBlock(i)
                mol_smi = Chem.MolToSmiles(i)
                try:
                    enumerator = MolStandardize.rdMolStandardize.TautomerEnumerator()
                    canonical_tautomer = enumerator.Canonicalize(mol_smi)
                    standard_mol_smi = Chem.MolToSmiles(canonical_tautomer)
                except:
                    standard_mol_smi = mol_smi
                
                molecule = None
                molecule = self.session.query(Molecule).filter(Molecule.smiles==standard_mol_smi).first()
                if molecule == None:
                    mol2 = Chem.MolFromSmiles(standard_mol_smi)
                    mol2 = Chem.AddHs(mol2)
                    mol_2d = Chem.MolToMolBlock(mol2)
                    molecule = Molecule(standard_mol_smi, mol_2d, 1)
                    self.session.add(molecule)
                    print('upload new molecule')
                    self.session.commit()
                    
                    check_molecule = self.session.query(Molecule).filter(Molecule.smiles==standard_mol_smi).first()
                    
                    if check_molecule !=None:
                        if check_molecule.conformer == []:
                            num_conformer = (0,)
                        else:
                            num_conformer = self.session.query(Conformer.num).filter_by(molecule_id=molecule.id).order_by(desc(Conformer.num)).first()
                        
                        if num_conformer[0] < 1001:
                            print(num_conformer)
                            new_conformer = Conformer('', 1001)
                            new_conformer.molecule = check_molecule
                        else:
                            new_conformer = Conformer('', num_conformer[0] + 1)
                            new_conformer.molecule = check_molecule
                    
                        self.session.add(new_conformer)
                        print('upload new conformer')
                        self.session.commit()
                        
                        job = self.session.query(Job).filter_by(id=self.job_id).first()
                        binding_site = self.session.query(Binding_site).filter_by(id=self.bsite_id).first()
                        result_conformer = Result_conformer(rmsd=0.0, mol_3d=mol_3d_block, score=glide_score)
                        result_conformer.hit = new_conformer
                        result_conformer.job_id = job.id
                        result_conformer.binding_site = binding_site
                    
                        self.session.add(result_conformer)
                        print('upload glide result_conf')
                        self.session.commit()
                        
                        save_dir = './media'
                   
                        
                        mol_img = Chem.MolFromSmiles(standard_mol_smi)
                        
                        image_path = os.path.join(save_dir, str(molecule.id) + '.png')
                        if not os.path.exists(image_path):
                            # print(image_path)
                            Draw.MolToFile(mol_img, image_path)
                            
                        bucket_path = "%s/%s/%s" % ('ncap_media', self.job_id, str(check_molecule.id) + '.png')
                        print(bucket_path)
                        blob = bucket.blob(bucket_path)
                        blob.upload_from_filename(image_path)
                        
                else:
                    new_conformer = Conformer(' ', 1001)
                    new_conformer.molecule = molecule
                    self.session.add(new_conformer)
                    print('upload new conformer')
                    self.session.commit()

                    job = self.session.query(Job).filter_by(id=self.job_id).first()
                    binding_site = self.session.query(Binding_site).filter_by(id=self.bsite_id).first()
                    result_conformer = Result_conformer(rmsd=0.0, mol_3d=mol_3d_block, score=glide_score)
                    result_conformer.hit = new_conformer
                    result_conformer.job_id = job.id
                    result_conformer.binding_site = binding_site
                
                    self.session.add(result_conformer)
                    print('upload glide result_conf')
                    self.session.commit()

                    save_dir = './media'
                  
                    
                    mol_img = Chem.MolFromSmiles(standard_mol_smi)
                    
                    image_path = os.path.join(save_dir, str(molecule.id) + '.png')
                    if not os.path.exists(image_path):
                        # print(image_path)
                        Draw.MolToFile(mol_img, image_path)
                        
                    bucket_path = "%s/%s/%s" % ('ncap_media', self.job_id, str(molecule.id) + '.png')
                    print(bucket_path)
                    blob = bucket.blob(bucket_path)
                    blob.upload_from_filename(image_path)
        
        
        
        
        time.sleep(10)
        
    
        
        file_path = f"./{self.file_name}_prep.maegz"
        bucket_path = "%s/%s" % ('PDB_prep', f"{self.file_name}_prep.maegz")
        #bucket_path = "PDB_prep/"
        print(bucket_path)
        blob = bucket.blob(bucket_path)
        blob.upload_from_filename(file_path)
    
    
        #====remove tmp make file====#
        os.remove(save_pdb)
        #os.remove(save_pdbp)
        os.remove(tmp_addH_file)
        os.remove(tmp_addH_rech)
