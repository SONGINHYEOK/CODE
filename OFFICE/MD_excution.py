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

class Molecule_dynamics:
    save_dir = '/home/NCAP_test/media'
    KST = 'Asia/Seoul'
    read_path = "PDB/"
    #read_path2 = "PDB_prep/"
    save_path = "./"
    
    def __init__(self, session, file_name, job_id):#, bsite_id):#, , job_id):
        self.session = session
        self.file_name = file_name
        self.job_id = job_id
        #self.bsite_id = bsite_id
        
    def finish_job(self):
        # self.session.rollback() # temp process - query feature included in session
        KST = pytz.timezone(self.KST)
        end_dt = datetime.now(KST)
        job = self.session.query(Job).filter_by(id=self.job_id).first()
        job.end_date = end_dt
        #self.session.commit()
        self.session.close()
    
    def act_md(self):
        '''
        cfg ="""
            annealing = false
            backend = {
            }
            bigger_rclone = false
            box = ?
            checkpt = {
            first = 0.0
            interval = 240.06
            name = "$JOBNAME.cpt"
            write_last_step = true
            }
            cpu = 1
            cutoff_radius = 9.0
            elapsed_time = 0.0
            energy_group = false
            eneseq = {
            first = 0.0
            interval = 1.2
            name = "$JOBNAME$[_replica$REPLICA$].ene"
            }
            ensemble = {
            barostat = {
                tau = 2.0
            }
            class = NPT
            method = MTK
            thermostat = {
                tau = 1.0
            }
            }
            glue = solute
            maeff_output = {
            center_atoms = solute
            first = 0.0
            interval = 120.0
            name = "$JOBNAME$[_replica$REPLICA$]-out.cms"
            periodicfix = true
            trjdir = "$JOBNAME$[_replica$REPLICA$]_trj"
            }
            meta = false
            meta_file = ?
            pressure = [1.01325 isotropic ]
            randomize_velocity = {
            first = 0.0
            interval = inf
            seed = 2007
            temperature = "@*.temperature"
            }
            restrain = none
            restraints = {
            existing = ignore
            new = []
            }
            simbox = {
            first = 0.0
            interval = 1.2
            name = "$JOBNAME$[_replica$REPLICA$]_simbox.dat"
            }
            surface_tension = 0.0
            taper = false
            temperature = [
            [300.0 0 ]
            ]
            time = 1200.0
            timestep = [0.002 0.002 0.006 ]
            trajectory = {
            center = []
            first = 0.0
            format = dtr
            frames_per_file = 250
            interval = 1.2
            name = "$JOBNAME$[_replica$REPLICA$]_trj"
            periodicfix = true
            write_last_vel = false
            write_velocity = false
            }
            """
        with open("MD_templete.cfg", 'w') as cf:
            cf.write(cfg)
           
        
        mjs = """
            # Desmond standard NPT relaxation protocol
            # All times are in the unit of ps.
            # Energy is in the unit of kcal/mol.
            task {
            task = "desmond:auto"
            set_family = {
                desmond = {
                    checkpt.write_last_step = no
                }
            }
            }

            simulate {
            title       = "Brownian Dynamics NVT, T = 10 K, small timesteps, and restraints on solute heavy atoms, 100ps"
            annealing   = off
            time        = 100
            timestep    = [0.001 0.001 0.003 ]
            temperature = 10.0
            ensemble = {
                class = "NVT"
                method = "Brownie"
                brownie = {
                    delta_max = 0.1
                }
            }
            restraints.new = [{
                name = posre_harm
                atoms = solute_heavy_atom
                force_constants = 50.0
            }]
            }

            simulate {
            title       = "NVT, T = 10 K, small timesteps, and restraints on solute heavy atoms, 12ps"
            annealing   = off
            time        = 12
            timestep    = [0.001 0.001 0.003]
            temperature = 10.0
            restraints.new = [{
                name = posre_harm
                atoms = solute_heavy_atom
                force_constants = 50.0
            }]
            ensemble = {
                class  = NVT
                method = Langevin
                thermostat.tau = 0.1
            }

            randomize_velocity.interval = 1.0
            eneseq.interval             = 0.3
            trajectory.center           = []
            }

            simulate {
            title       = "NPT, T = 10 K, and restraints on solute heavy atoms, 12ps"
            annealing   = off
            time        = 12
            temperature = 10.0
            restraints.existing = retain
            ensemble    = {
                class  = NPT
                method = Langevin
                thermostat.tau = 0.1
                barostat  .tau = 50.0
            }

            randomize_velocity.interval = 1.0
            eneseq.interval             = 0.3
            trajectory.center           = []
            }

            simulate {
            title       = "NPT and restraints on solute heavy atoms, 12ps"
            effect_if   = [["@*.*.annealing"] 'annealing = off temperature = "@*.*.temperature[0][0]"']
            time        = 12
            restraints.existing = retain
            ensemble    = {
                class  = NPT
                method = Langevin
                thermostat.tau = 0.1
                barostat  .tau = 50.0
            }

            randomize_velocity.interval = 1.0
            eneseq.interval             = 0.3
            trajectory.center           = []
            }

            simulate {
            title       = "NPT and no restraints, 24ps"
            effect_if   = [["@*.*.annealing"] 'annealing = off temperature = "@*.*.temperature[0][0]"']
            time        = 24
            ensemble    = {
                class  = NPT
                method = Langevin
                thermostat.tau = 0.1
                barostat  .tau = 2.0
            }

            eneseq.interval   = 0.3
            trajectory.center = solute
            }

            simulate {
            cfg_file = "./MD_templete.cfg"
            jobname  = "$MAINJOBNAME"
            dir      = "."
            compress = ""
            }

            pl_analysis {
                ligand_asl = ""
                protein_asl = ""
            }

            """
        with open(f"{self.file_name}_md.msj", 'w') as mj:
            mj.write(mjs)
        
        
    
        start_md(self.file_name)
        print('start MD')
        while os.path.exists('/home/NCAP_test/'+f"{self.file_name}_MD.eaf") !=True:
            time.sleep(10)
            if os.path.exists('/home/NCAP_test/'+f"{self.file_name}_MD.eaf") ==True:
                break
        print('complete MD')
        
   
        print('start align')
        align_cms(self.file_name)
        print('end align')
    
        print('start extract dat')
        ext_DAT(self.file_name)
        
        while os.path.exists('./PL_RMSD.dat') !=True:
            time.sleep(10)
            if os.path.exists('./PL_RMSD.dat') ==True:
                break
        print('end extract dat')
     
        
        
        print('start extract PDB,XTC')
        os.system(f"/apps/schrodinger2022-1/run /home/NCAP_test/trj_no_virt.py --no-h -o {self.file_name}_MD_result align_{self.file_name}-out.cms align_{self.file_name}_trj")
        print('end extract PDB,XTC')
        
        lig_asl=ext_asl(f"align_{self.file_name}-out.cms")
        
        print("start MMGBSA")
        start_mmgbsa(lig_asl,f"mmgbsa_{self.file_name}" ,f"align_{self.file_name}-out.cms" )
        while os.path.exists('/home/NCAP_test/'+f"mmgbsa_{self.file_name}-prime-out.csv") !=True:
            time.sleep(10)
            if os.path.exists('/home/NCAP_test/'+f"mmgbsa_{self.file_name}-prime-out.csv") ==True:
                break
        print("end MMGBSA")

        '''
       
        gce_instance = 'inhyeok.song@nt-ddp-ncap-web-server'
        gce_project = 'nettargets-chemdb-v1'
        gce_zone = 'asia-northeast3-b'
        gce_query_path = f'/home/inhyeok.song/Ntrophy/media/MD_result/{self.job_id}'
        
        gce_output_path = f'/home/NCAP_test/{self.file_name}_MD_result.pdb'
        gce_output_path1 = f'/home/NCAP_test/{self.file_name}_MD_result.xtc'
        gce_output_path2 = f'/home/NCAP_test/mmgbsa_{self.file_name}-prime-out.csv'
        gce_output_path3 = f'/home/NCAP_test/{self.file_name}_MD.ene'
        gce_output_path4 = '/home/NCAP_test/PL_RMSD.dat'
        gce_output_path5 = '/home/NCAP_test/PL-Contacts_HBond.dat'
        gce_output_path6 = '/home/NCAP_test/PL-Contacts_Hydrophobic.dat'
        gce_output_path7 = '/home/NCAP_test/PL-Contacts_Ionic.dat'
        gce_output_path8 = '/home/NCAP_test/PL-Contacts_WaterBridge.dat'
        gce_output_path9 = '/home/NCAP_test/P_RMSF.dat'
        #gce_output_path10 = f'/home/NCAP_test/mmgbsa_{self.file_name}-prime-out.csv'
        total_list = [gce_output_path, gce_output_path1, gce_output_path2, gce_output_path3, gce_output_path4, gce_output_path5, gce_output_path6, gce_output_path7,gce_output_path8, gce_output_path9]#, gce_output_path10]
        
        for total in total_list:    
            cmd_scp = "gcloud compute scp %s %s:%s --zone %s --project %s" % (total, gce_instance, gce_query_path, gce_zone, gce_project)
            process = subprocess.Popen(cmd_scp.split(), stdout=subprocess.PIPE).wait(timeout=None)
        
