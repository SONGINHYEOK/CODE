import os
from os import listdir
from os.path import isfile, join
import time
import pandas as pd
import numpy as np
import multiprocessing
from schprocess import *
from datetime import date
from multiprocessing import Process, current_process
from dbmodel import Compound, Molecule, Compound_Molecule
from dbmodel import Conformer, ConformerFingerprint
from dbmodel import Pharmacophore, PharmacophoreType, PharmacophoreDistance, DistanceType
from dbmodel import Protein_structure, Ligand, Binding_site, Site_residue, Amino_acid
from dbmodel import Molecule_fragment, Fragment, Fragment_edge, Scaffold, Reagent, Reagent_compound
from dbmodel import Job, Result_conformer
from featureutil import make_pair_distance, make_distance_fingerprint, parse_feature_sch
from extractprotein import *
from standardizeutil import *
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

class ScaffoldUploader:
    def __init__(self, session, id_list):
        self.session = session
        self.id_list = id_list

    def upload_scaffold(self):
        proc_name = current_process().name
        self.process_scaffold()
        self.session.close()

    def process_scaffold(self):
        # alr_mol = self.session.query(Molecule_fragment).order_by(desc(Molecule_fragment.molecule_id)).first()
        # print(alr_mol.molecule_id)
        # molecules = self.session.query(Molecule).filter(Molecule.id > alr_mol.molecule_id).order_by(Molecule.id).all()

        n = 0
        if self.id_list == None:
            #molecules = self.session.query(Molecule).join(Conformer, Molecule.id == Conformer.molecule_id).join(Result_conformer, Result_conformer.hit_id == Conformer.id).filter(Result_conformer.job_id == 6).all()
            df = pd.read_csv("/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/NT006/92_final_molecule.csv")
            id_list = df['id'].tolist()
            #id_list =["16430513"]
            molecules = self.session.query(Molecule).filter(Molecule.id.in_(id_list)).all()
            #molecules = self.session.query(Molecule).order_by(Molecule.id).all()
            #molecules = self.session.query(Molecule).order_by(Molecule.id).limit().all()
        else:
            molecules = self.session.query(Molecule).filter(Molecule.id.in_(self.id_list)).all()

        for molecule in molecules:
            n += 1
            start_time = time.time()
            smi = molecule.smiles
            if len(smi) <415:
                #already_mol = None
                #already_mol = self.session.query(Molecule_fragment).filter_by(molecule_id = molecule.id).first()
                #if already_mol != None:
                #    continue
                
                molecule_fragment = molecule.molecule_fragment
 
                
                fragment_boolean = False
                scaffold_boolean = False
                
                            
                if len(molecule_fragment) > 0:
                    fragment_boolean = True
                    for mol_frag in molecule_fragment:
                        if mol_frag.scaffold != None:
                            scaffold_boolean = True
                            break
                
                if fragment_boolean == False and scaffold_boolean==False:
                    print('not at all')
                    mol = Chem.MolFromSmiles(smi)

                    print(molecule.id, smi)
                    params = rdScaffoldNetwork.ScaffoldNetworkParams()
                    params.includeScaffoldsWithAttachments=False
                    try:
                        netwks = rdScaffoldNetwork.CreateScaffoldNetwork([mol], params)
                    except:
                        continue
                    
                    node_list = []
                    edge_list = []
                    frag_list = {}
                    
                    for edge in netwks.edges:
                        if edge.beginIdx != 0:
                            edge_list.append([edge.beginIdx - 1, edge.endIdx - 1 ])

                    
                    for i, node in enumerate(netwks.nodes):
                        if i != 0:
                            node_list.append([node, netwks.edges[i-1].type])
                   
                    
                    for i, node in enumerate(node_list):
                        try:
                            (molSmiles, neutralised) = NeutraliseCharges(node[0])
                            standard_smi = MolStandardize.canonicalize_tautomer_smiles(molSmiles)
                        except:
                            break

                        frag = None
                        frag = self.session.query(Fragment).filter_by(smiles=str(standard_smi)).first()
                        if not frag:
                            if str(node[1]) == 'Generic': 
                                frag = Fragment(smiles=standard_smi, generic=True)
                            else:
                                frag = Fragment(smiles=standard_smi, generic=False)
                            self.session.add(frag)
                            print("add fragment")
                            self.session.commit()
                            frag_list[i] = [frag, False]

                            # mirroring
                            #conn_info = "host='35.220.207.27' dbname='chemdb' user='postgres' password='1234' port='5433'"
                            #conn = psycopg2.connect(conn_info)
                            #cur = conn.cursor()

                            #command = "insert into mirror_fragment(id, canonical_smiles, generic) values('" + str(frag.id) + "','" + str(frag.smiles) + "','" + str(frag.generic) + "');"
                            #cur.execute(command)
                            #conn.commit()
                            #cur.close()
                            #conn.close()

                        else :
                            frag_list[i] = [frag , True]

                        # add molecule_fragment
                        molecule_fragment = None
                        molecule_fragment = self.session.query(Molecule_fragment).filter(Molecule_fragment.fragment_id == frag.id, Molecule_fragment.molecule_id == molecule.id).first()
                        if not molecule_fragment:
                            molecule_fragment = Molecule_fragment()
                            molecule_fragment.fragment = frag
                            molecule_fragment.molecule = molecule
                            self.session.add(molecule_fragment)
                            print("add molecule_fragment")
                            self.session.commit()

                        # add scaffold
                        if str(node[1]) == 'Initialize':
                            scaffold = None
                            scaffold = self.session.query(Scaffold).filter(Scaffold.mol_frag_id == molecule_fragment.id).first()
                            if not scaffold:
                                scaffold = Scaffold()
                                scaffold.molecule_fragment = molecule_fragment
                                self.session.add(scaffold)           
                                print("add scaffold")   
                                self.session.commit()

                        if frag_list == {} :
                            continue
                        else:
                            for i, edge in enumerate(edge_list):
                                print(frag_list, edge)
                                
                                fragment_st = frag_list[edge[0]]
                                fragment_end = frag_list[edge[1]]
                               
                                if fragment_st[1] == True and fragment_end[1] == True:
                                    print("already fragment_edge")
                                    continue
                                else:
                                    fragment_edge = Fragment_edge()
                                    fragment_edge.fragment_st = fragment_st[0]
                                    fragment_edge.fragment_end = fragment_end[0]
                                    self.session.add(fragment_edge)
                                    self.session.commit()
                                    print("add fragment edge")                        
                elif fragment_boolean == True and scaffold_boolean == False:
                    print(molecule.id, smi)
                    mol = Chem.MolFromSmiles(smi)

                    params = rdScaffoldNetwork.ScaffoldNetworkParams()
                    params.includeScaffoldsWithAttachments=False
                    
              
                    try:
                        netwks = rdScaffoldNetwork.CreateScaffoldNetwork([mol], params)
                    except:
                        continue
                    
                    node_list = []
                     
                    for i, node in enumerate(netwks.nodes):
                        if i != 0:
                            node_list.append([node, netwks.edges[i-1].type])                 
                        if str(netwks.edges[i-1].type) == 'Initialize':
                            molecule_fragment = molecule.molecule_fragment
                            
                            scaffold = None
                            scaffold = self.session.query(Scaffold).filter(Scaffold.mol_frag_id == molecule_fragment.id).first()
                            if not scaffold:
                                scaffold = Scaffold()
                                scaffold.molecule_fragment = molecule_fragment
                                self.session.add(scaffold)           
                                print("add scaffold")   
                                self.session.commit()

                            for i, edge in enumerate(edge_list):
                                fragment_st = frag_list[edge[0]]
                                fragment_end = frag_list[edge[1]]

                                if fragment_st[1] == True and fragment_end[1] == True:
                                    print("already fragment_edge")
                                    continue
                                else:
                                    fragment_edge = Fragment_edge()
                                    fragment_edge.fragment_st = fragment_st[0]
                                    fragment_edge.fragment_end = fragment_end[0]
                                    self.session.add(fragment_edge)
                                    self.session.commit()
                                    print("add fragment edge")
                            
                    if 'Initialize' not in node_list:
                        frag = None
                        frag = self.session.query(Fragment).filter_by(smiles=str(node_list[0][0])).first()
                        
                        if frag == None:
                            frag = Fragment(smiles=str(node_list[0][0]), generic=False)
                            self.session.add(frag)
                            print("not intialize add fragment")
                            self.session.commit()
            
                            
                            molecule_fragment = None
                            molecule_fragment = self.session.query(Molecule_fragment).filter(Molecule_fragment.fragment_id == frag.id, Molecule_fragment.molecule_id == molecule.id).first()
                            if not molecule_fragment:
                                molecule_fragment = Molecule_fragment()
                                molecule_fragment.fragment = frag
                                molecule_fragment.molecule = molecule
                                self.session.add(molecule_fragment)
                                print("not intialize add molecule_fragment")
                                self.session.commit()

                            scaffold = None
                            scaffold = self.session.query(Scaffold).filter(Scaffold.mol_frag_id == molecule_fragment.id).first()
                            if not scaffold:
                                scaffold = Scaffold()
                                scaffold.molecule_fragment = molecule_fragment
                                self.session.add(scaffold)           
                                print("not intialize add scaffold")   
                                self.session.commit()
                        else:
                                                        
                            molecule_fragment = None
                            molecule_fragment = self.session.query(Molecule_fragment).filter(Molecule_fragment.fragment_id == frag.id, Molecule_fragment.molecule_id == molecule.id).first()
                            if not molecule_fragment:
                                molecule_fragment = Molecule_fragment()
                                molecule_fragment.fragment = frag
                                molecule_fragment.molecule = molecule
                                self.session.add(molecule_fragment)
                                print("not intialize add molecule_fragment")
                                self.session.commit()

                            scaffold = None
                            scaffold = self.session.query(Scaffold).filter(Scaffold.mol_frag_id == molecule_fragment.id).first()
                            if not scaffold:
                                scaffold = Scaffold()
                                scaffold.molecule_fragment = molecule_fragment
                                self.session.add(scaffold)           
                                print("not intialize add scaffold")   
                                self.session.commit()
                                
                else:
                    continue
       
class Scaffold_Heavy_Uploader:
    def __init__(self, SessionFactory):
        self.SessionFactory = SessionFactory
        self.session = self.SessionFactory()
        self.proc_size = 20

    def upload_scaffold_multi(self):
        #mol_id_list = self.session.query(Molecule.id).all()
        mol_id_list = self.session.query(Molecule.id).order_by(Molecule.id).limit(1000000).all()
        mol_id_list = [ mol[0] for mol in mol_id_list ]

        id_list = []
        for i in range(0, len(mol_id_list), 1000):
            id_list.append(mol_id_list[i : i+1000])

        self.item_to_proc(id_list, 'scaffold')

    def item_to_proc(self, num_list, type):
        que = multiprocessing.JoinableQueue()
        procs = []

        for i in range(self.proc_size):
            proc = Process(target=self.provide_queue, args=(que, type))
            procs.append(proc)
            time.sleep(0.1)
            procs[-1].daemon = True
            procs[-1].start()
        
        for num in num_list:
            que.put(num)
        que.join()
        
        for proc in procs:
            que.put(None)
        que.join()

        for proc in procs:
            proc.join
    
    def provide_queue(self, que, type):
        for item in iter(que.get, None):
            self.upload(item, type)
            que.task_done()
        que.task_done()

    def upload(self, item, type):
        session = self.SessionFactory()
        uploader = ScaffoldUploader(session, item)
        uploader.upload_scaffold()