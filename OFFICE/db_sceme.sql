---
PostgreSQL

-- ---
-- Globals
-- ---

-- SET SQL_MODE="NO_AUTO_VALUE_ON_ZERO";
-- SET FOREIGN_KEY_CHECKS=0;

-- ---
-- Table 'compound'
-- 
-- ---

DROP TABLE IF EXISTS compound;
CREATE TABLE compound (
  id bigserial primary key,
  source_id INTEGER NOT NULL DEFAULT NULL,
  original_smiles VARCHAR(10000) NOT NULL,
  original_id VARCHAR(100) NOT NULL,
  created_at date NOT NULL,
  updated_at date NOT NULL
);

-- ---
-- Table 'compound_molecule'
-- 
-- ---

DROP TABLE IF EXISTS compound_molecule;
CREATE TABLE compound_molecule (
  id bigserial primary key,
  compound_id BIGINT NOT NULL,
  molecule_id BIGINT NOT NULL
);

-- ---
-- Table 'molecule'
-- 
-- ---

DROP TABLE IF EXISTS molecule;
CREATE TABLE molecule (
  id bigserial primary key,
  canonical_smiles VARCHAR(10000) NOT NULL,
  inchi VARCHAR(10000) DEFAULT NULL,
  inchi_key VARCHAR(50) DEFAULT NULL,
  mol_2d VARCHAR(1000000) NOT NULL,
  path VARCHAR(100)
);

-- ---
-- Table 'source'
-- 
-- ---

DROP TABLE IF EXISTS source;
CREATE TABLE source (
  id bigserial primary key,
  name VARCHAR(50) NOT NULL
);

-- ---
-- Table 'vendor'
-- 
-- ---

DROP TABLE IF EXISTS vendor;
CREATE TABLE vendor (
  id bigserial primary key,
  name VARCHAR(100) NOT NULL
);

-- ---
-- Table 'reagent'
-- 
-- ---

DROP TABLE IF EXISTS reagent;
CREATE TABLE reagent (
  id bigserial primary key,
  vendor_id INTEGER NOT NULL,
  code VARCHAR(100) NOT NULL
);

-- ---
-- Table 'reagent_compound'
-- 
-- ---

DROP TABLE IF EXISTS reagent_compound;
CREATE TABLE reagent_compound (
  id bigserial primary key,
  compound_id BIGINT NOT NULL,
  reagent_id BIGINT NOT NULL
);

-- ---
-- Table 'conformer'
-- 
-- ---

DROP TABLE IF EXISTS conformer;
CREATE TABLE conformer (
  id bigserial primary key,
  molecule_id BIGINT NOT NULL,
  conformer_number INTEGER NOT NULL
);

-- ---
-- Table 'usr'
-- 
-- ---

DROP TABLE IF EXISTS member;
CREATE TABLE member (
  id serial primary key,
  name VARCHAR(50) NOT NULL,
  password VARCHAR(50) NOT NULL,
  user_id VARCHAR(20) NULL,
  role_id VARCHAR(30) NULL,
  use_yn VARCHAR(1) NULL DEFAULT 'N',
  creat_user VARCHAR(30) NULL,
  creat_dt timestamptz(0) NULL,
  updt_user VARCHAR(30) NULL,
  updt_dt timestamptz(0) NULL
);

-- -- ---
-- -- Table 'pharmacophore'
-- -- 
-- -- ---

-- DROP TABLE IF EXISTS pharmacophore;
-- CREATE TABLE pharmacophore (
--   id bigserial primary key,
--   conformer_id BIGINT NOT NULL,
--   type_id INTEGER NOT NULL,
--   head_x FLOAT DEFAULT NULL,
--   head_y FLOAT DEFAULT NULL,
--   head_z FLOAT DEFAULT NULL,
--   head_tolerance FLOAT DEFAULT NULL,
--   tail_x FLOAT NOT NULL,
--   tail_y FLOAT NOT NULL,
--   tail_z FLOAT NOT NULL,
--   tail_tolerance FLOAT NOT NULL
-- );

-- ---
-- Table 'pharmacophore_type'
-- 
-- ---

DROP TABLE IF EXISTS pharmacophore_type;
CREATE TABLE pharmacophore_type (
  id serial primary key,
  name VARCHAR(50) NOT NULL
);

-- ---
-- Table 'pharmacophore_distance'
--
-- ---

--DROP TABLE IF EXISTS pharmacophore_distance;
--CREATE TABLE pharmacophore_distance (
--  id bigserial primary key,
--  pharmacophore1_id BIGINT NOT NULL,
--  pharmacophore2_id BIGINT NOT NULL,
--  type_id INTEGER NOT NULL,
--  distance FLOAT NOT NULL
--);

-- ---
-- Table 'distance_type'
--
-- ---

--DROP TABLE IF EXISTS distance_type;
--CREATE TABLE distance_type (
--  id serial primary key,
--  name VARCHAR(30) NOT NULL
--);

-- ---
-- Table 'conformer_fingerprint'
-- 
-- ---

-- DROP TABLE IF EXISTS conformer_fingerprint;
-- CREATE TABLE conformer_fingerprint (
--   id BIGINT primary key,
--   fp BIGINT NOT NULL,
--   aa BIGINT NOT NULL,
--   ad BIGINT NOT NULL,
--   ah BIGINT NOT NULL,
--   an BIGINT NOT NULL,
--   ap BIGINT NOT NULL,
--   ar BIGINT NOT NULL,
--   dd BIGINT NOT NULL,
--   dh BIGINT NOT NULL,
--   dn BIGINT NOT NULL,
--   dp BIGINT NOT NULL,
--   dr BIGINT NOT NULL,
--   hh BIGINT NOT NULL,
--   hn BIGINT NOT NULL,
--   hp BIGINT NOT NULL,
--   hr BIGINT NOT NULL,
--   nn BIGINT NOT NULL,
--   np BIGINT NOT NULL,
--   nr BIGINT NOT NULL,
--   pp BIGINT NOT NULL,
--   pr BIGINT NOT NULL,
--   rr BIGINT NOT NULL
-- );

-- ---
-- Table 'ligand'
-- 
-- ---

DROP TABLE IF EXISTS ligand;

CREATE TABLE ligand (
  id serial primary key,
  ligand_id VARCHAR(10) NOT NULL, 
  molecule_id BIGINT NOT NULL,
  mol_3d VARCHAR(1000000) NOT NULL
);

-- ---
-- Table 'binding_site'
-- 
-- ---

DROP TABLE IF EXISTS binding_site;
CREATE TABLE binding_site (
  id serial primary key,
  protein_str_id INTEGER NOT NULL,
  ligand_id INTEGER NULL
);

-- ---
-- Table 'site_residue'
-- 
-- ---

DROP TABLE IF EXISTS site_residue;
CREATE TABLE site_residue (
  id serial primary key,
  chain CHAR(1) NOT NULL,
  residue_number INTEGER NOT NULL,
  residue_type_id INTEGER NOT NULL,
  binding_site_id INTEGER NOT NULL
);

-- ---
-- Table 'protein_structure'
-- 
-- ---

DROP TABLE IF EXISTS protein_structure;

CREATE TABLE protein_structure (
  id serial primary key,
  pdb_id CHAR(50) NULL DEFAULT NULL,
  source_id INTEGER NOT NULL,
  name VARCHAR(500) NOT NULL,
  path VARCHAR(200) NOT NULL
);

-- ---
-- Table 'amino_acid'
-- 
-- ---

DROP TABLE IF EXISTS amino_acid;

CREATE TABLE amino_acid (
  id serial primary key,
  synonym CHAR(30) NOT NULL,
  abbreviation CHAR(3) NOT NULL,
  name VARCHAR(50) NOT NULL
);

-- ---
-- Table 'fragment'
--
-- ---

DROP TABLE IF EXISTS fragment;
CREATE TABLE fragment (
  id bigserial primary key,
  canonical_smiles VARCHAR(10000) NOT NULL,
  generic BOOLEAN NOT NULL DEFAULT '0'
);

-- ---
-- Table 'molecule_fragment'
-- 
-- ---

DROP TABLE IF EXISTS molecule_fragment;
CREATE TABLE molecule_fragment (
  id bigserial primary key,
  molecule_id BIGINT NOT NULL,
  fragment_id BIGINT NOT NULL
);

-- ---
-- Table 'fragment_edge'
-- 
-- ---

DROP TABLE IF EXISTS fragment_edge;
CREATE TABLE fragment_edge (
  id bigserial primary key,
  start BIGINT NOT NULL,
  "end" BIGINT NOT NULL
);

-- ---
-- Table 'scaffold'
-- 
-- ---

DROP TABLE IF EXISTS scaffold;
CREATE TABLE scaffold (
  id bigserial primary key,
  mol_frag_id BIGINT NOT NULL
);

-- ---
-- Table 'virtual'
-- 
-- ---

DROP TABLE IF EXISTS virtual;
CREATE TABLE virtual (
  id BIGINT NOT NULL PRIMARY KEY
);



-- ---
-- Table 'cluster'
-- 
-- ---

DROP TABLE IF EXISTS cluster;
CREATE TABLE cluster (
  id serial primary key,
  number INTEGER NOT NULL,
  fragment_id BIGINT NULL DEFAULT NULL
);


-- ---
-- Table 'custom'
-- 
-- ---

--DROP TABLE IF EXISTS custom;
--CREATE TABLE custom (
--  id INTEGER PRIMARY KEY
--);

-- ---
-- Table 'protein_target'
--
-- ---
DROP TABLE IF EXISTS protein_target;
CREATE TABLE protein_target (
 id serial primary key,
 code VARCHAR NOT NULL,
 function VARCHAR NOT NULL,
 disruption_phenotype VARCHAR NULL DEFAULT NULL,
 tissue_specificity VARCHAR NULL DEFAULT NULL,
 sequence VARCHAR NOT NULL,
 gene_id INTEGER NULL DEFAULT NULL,
 chembl_id VARCHAR NULL DEFAULT NULL
);

-- ---
-- Table 'uniprot_code'
-- 
-- ---

DROP TABLE IF EXISTS uniprot_code;
CREATE TABLE uniprot_code (
  id serial primary key,
  accession VARCHAR(50) NOT NULL
);

-- ---
-- Table 'protein_name'
-- 
-- ---

DROP TABLE IF EXISTS protein_name;
CREATE TABLE protein_name (
  id serial primary key,
  protein_id INTEGER NOT NULL,
  name VARCHAR(500) NOT NULL
);

-- ---
-- Table 'protein_pdb'
-- 
-- ---

DROP TABLE IF EXISTS protein_pdb;
CREATE TABLE protein_pdb (
  id serial primary key,
  protein_id INTEGER NOT NULL,
  protein_str_id INTEGER NOT NULL,
  method INTEGER NOT NULL,
  resolution FLOAT NULL DEFAULT NULL
);

-- ---
-- Table 'gene_name'
-- 
-- ---

DROP TABLE IF EXISTS gene_name;
CREATE TABLE gene_name (
  id serial primary key,
  protein_id INTEGER NOT NULL,
  name VARCHAR(50) NOT NULL
);

-- ---
-- Table 'organism'
-- 
-- ---

DROP TABLE IF EXISTS organism;
CREATE TABLE organism (
  id serial primary key,
  taxonomy_id INTEGER NOT NULL,
  scientific_name VARCHAR(500) NOT NULL,
  common_name VARCHAR(100) NULL DEFAULT NULL,
  synonym VARCHAR(100) NULL DEFAULT NULL
);

-- ---
-- Table 'protein_organism'
-- 
-- ---

DROP TABLE IF EXISTS protein_organism;
CREATE TABLE protein_organism (
  id serial primary key,
  protein_id INTEGER NOT NULL,
  organism_id INTEGER NOT NULL
);

-- ---
-- Table 'disease'
-- 
-- ---

DROP TABLE IF EXISTS disease;
CREATE TABLE disease (
  id serial primary key,
  accession VARCHAR(20) NOT NULL,
  name VARCHAR(500) NOT NULL,
  acronym VARCHAR(20) NOT NULL,
  definition VARCHAR(10000) NOT NULL
);

-- ---
-- Table 'protein_disease'
-- 
-- ---

DROP TABLE IF EXISTS protein_disease;
CREATE TABLE protein_disease (
  id serial primary key,
  protein_id INTEGER NOT NULL,
  disease_id INTEGER NOT NULL
);

-- ---
-- Table 'keyword'
-- 
-- ---

DROP TABLE IF EXISTS keyword;
CREATE TABLE keyword (
  id serial primary key,
  accession VARCHAR(50) NULL,
  name VARCHAR(100) NOT NULL,
  definition VARCHAR(10000) NOT NULL,
  go VARCHAR(500) NULL DEFAULT NULL,
  category INTEGER NULL DEFAULT NULL
);

-- ---
-- Table 'protein_keyword'
-- 
-- ---

DROP TABLE IF EXISTS protein_keyword;
CREATE TABLE protein_keyword (
  id serial primary key,
  protein_id INTEGER NOT NULL,
  keyword_id INTEGER NOT NULL
);

-- ---
-- Table 'ppi'
-- 
-- ---

DROP TABLE IF EXISTS ppi;
CREATE TABLE ppi (
  id serial primary key,
  start INTEGER NOT NULL,
  "end" INTEGER NOT NULL
);

-- ---
-- Table 'ref_code'
-- 
-- ---

DROP TABLE IF EXISTS ref_code;
CREATE TABLE ref_code (
  id serial primary key,
  code VARCHAR(10) NOT NULL,
  source_id INTEGER NOT NULL
);

-- ---
-- Table 'ref_source'
-- 
-- ---

DROP TABLE IF EXISTS ref_source;
CREATE TABLE ref_source (
  id serial primary key,
  name VARCHAR(10) NOT NULL
);

-- ---
-- Table 'disease_reference'
-- 
-- ---

DROP TABLE IF EXISTS disease_reference;
CREATE TABLE disease_reference (
  id serial primary key,
  disease_id INTEGER NOT NULL,
  reference_id INTEGER NOT NULL
);

-- ---
-- Table 'disease_synonym'
-- 
-- ---

DROP TABLE IF EXISTS disease_synonym;
CREATE TABLE disease_synonym (
  id serial primary key,
  name VARCHAR(500) NOT NULL DEFAULT 'NULL',
  disease_id INTEGER NOT NULL
);

-- ---
-- Table 'disease_keyword'
-- 
-- ---

DROP TABLE IF EXISTS disease_keyword;
CREATE TABLE disease_keyword (
  id serial primary key,
  disease_id INTEGER NOT NULL,
  keyword_id INTEGER NOT NULL
);

-- ---
-- Table 'protein_uniprot'
-- 
-- ---

DROP TABLE IF EXISTS protein_uniprot;
CREATE TABLE protein_uniprot (
  id serial primary key, 
  uniprot_id INTEGER NOT NULL, 
  protein_id INTEGER NOT NULL
);

-- ---
-- Table 'pdb_chain'
-- 
-- ---

DROP TABLE IF EXISTS pdb_chain;
CREATE TABLE pdb_chain (
id serial primary key, 
chains VARCHAR(100) NOT NULL, 
start INTEGER NOT NULL, 
"end" INTEGER NOT NULL, 
protein_pdb_id INTEGER NOT NULL
);

-- ---
-- Table 'pdb_method'
--
-- ---

DROP TABLE IF EXISTS pdb_method;
CREATE TABLE pdb_method (
 id serial PRIMARY KEY,
 name VARCHAR(10) NOT NULL
);

-- ---
-- Table 'project'
-- 
-- ---

DROP TABLE IF EXISTS project;
CREATE TABLE project (
  id serial primary key, 
  name VARCHAR(200) NOT NULL,
  description varchar(1000) NULL,
  ntce_bgnde date NULL,
  ntce_endde date NULL,
  creat_user varchar(30) NULL,
  creat_dt timestamp NULL,
  updt_user varchar(30) NULL,
  updt_dt timestamp NULL,
  status varchar(10) NULL
);

-- ---
-- Table 'project_member'
-- 
-- ---

DROP TABLE IF EXISTS project_member;
CREATE TABLE project_member(
  id serial primary key, 
  project_id INTEGER NOT NULL,
  member_id INTEGER NOT NULL,
  project_role_id VARCHAR(30) NULL
);


-- ---
-- Table 'batch'
-- 
-- ---

DROP TABLE IF EXISTS batch;
CREATE TABLE batch (
  id serial primary key, 
  project_id INTEGER NOT NULL,
  name VARCHAR(10) NOT NULL,
  member_id INTEGER NULL,
  target_id INTEGER NULL
);

-- ---
-- Table 'job'
-- 
-- ---

DROP TABLE IF EXISTS job;
CREATE TABLE job (
  id serial primary key,
  project_member_id INTEGER NOT NULL,
  batch_id INTEGER NOT NULL,
  type_id INTEGER NOT NULL,
  start_date TIMESTAMP WITH TIME ZONE NOT NULL,
  end_date TIMESTAMP WITH TIME ZONE NULL DEFAULT NULL,
  job_name VARCHAR(50) NULL
);

-- ---
-- Table 'job_type'
--
-- ---

DROP TABLE IF EXISTS job_type;
CREATE TABLE job_type(
  id serial primary key,
  name VARCHAR(100) NOT NULL
);

-- ---
-- Table 'job_edge'
-- 
-- ---

DROP TABLE IF EXISTS job_edge;
CREATE TABLE job_edge (
  id serial primary key, 
  start INTEGER NULL,
  "end" INTEGER NOT NULL,
  batch_id INTEGER
);

-- ---
-- Table 'result_conformer'
-- 
-- ---

DROP TABLE IF EXISTS result_conformer;
CREATE TABLE result_conformer (
  id serial primary key,
  job_id INTEGER NOT NULL,
  hit_id BIGINT NOT NULL,
  bsite_id BIGINT NOT NULL,
  rmsd FLOAT NOT NULL,
  mol_3d VARCHAR(1000000) NOT NULL,
  score FLOAT,
  previous_yn BOOLEAN DEFAULT false
);

-- ---
-- Table 'result_molecule'
-- 
-- ---

DROP TABLE IF EXISTS result_molecule;
CREATE TABLE result_molecule (
  id serial primary key,
  job_id INTEGER NOT NULL,
  cluster_id INTEGER NULL DEFAULT NULL,
  molecule_id BIGINT NOT NULL,
  previous_yn BOOLEAN DEFAULT false
);

-- ---
-- Table 'result_assay'
--
-- ---

DROP TABLE IF EXISTS result_assay;
CREATE TABLE result_assay (
  id serial primary key,
  job_id INTEGER NOT NULL,
  molecule_id INTEGER NOT NULL,
  value FLOAT NOT NULL,
  previous_yn BOOLEAN DEFAULT false
);

-- ---
-- Table 'lead'
-- 
-- ---

DROP TABLE IF EXISTS lead;
CREATE TABLE lead(
  id serial primary key,
  job_conformer_id INTEGER NOT NULL
  --result_assay_id INTEGER NOT NULL
);

-- ---
-- Table 'candidate'
-- 
-- ---

DROP TABLE IF EXISTS candidate;
CREATE TABLE candidate(
  id serial primary key,
  job_conformer_id INTEGER NOT NULL
  --result_assay_id INTEGER NOT NULL
);


-- ---
-- Table 'job_molecule'
-- 
-- ---

--DROP TABLE IF EXISTS job_molecule;
--CREATE TABLE job_molecule (
--id serial primary key,
--project_id INTEGER NOT NULL,
--type_id INTEGER NOT NULL,
--start_date TIMESTAMP WITH TIME ZONE NOT NULL,
--end_date TIMESTAMP WITH TIME ZONE NULL DEFAULT NULL
--);

-- ---
-- Table 'job_conformer'
-- 
-- ---

--DROP TABLE IF EXISTS job_conformer;
--CREATE TABLE job_conformer (
--id serial primary key,
--project_id INTEGER NOT NULL,
--protein VARCHAR(500) NOT NULL,
--start_date TIMESTAMP WITH TIME ZONE NOT NULL,
--end_date TIMESTAMP WITH TIME ZONE NULL DEFAULT NULL,
--type INTEGER NOT NULL,
--name VARCHAR(50)
--);

-- ---
-- Table 'job_molecule_type'
-- 
-- ---

--DROP TABLE IF EXISTS job_molecule_type;
--CREATE TABLE job_molecule_type (
--id serial primary key,
--name VARCHAR(100) NOT NULL
--);

-- ---
-- Table 'job_conformer_type'
-- 
-- ---

--DROP TABLE IF EXISTS job_conformer_type;
--CREATE TABLE job_conformer_type (
--id serial primary key,
--name VARCHAR(100) NOT NULL,
--);

-- ---
-- Table 'dataset'
-- 
-- ---

-- DROP TABLE IF EXISTS dataset;
-- CREATE TABLE dataset (
--   id serial primary key,
--   name VARCHAR(200) NOT NULL
-- );

-- ---
-- Table 'dataset_molecule'
-- 
-- ---

-- DROP TABLE IF EXISTS dataset_molecule;
--   CREATE TABLE dataset_molecule (
--   id serial primary key,
--   dataset_id INTEGER NOT NULL,
--   molecule_id BIGINT NOT NULL
-- );

-- ---
-- Table ‘nmap_model’
-- 
-- ---
DROP TABLE IF EXISTS nmap_model;
CREATE TABLE nmap_model (
  id serial primary key,
  disease_id INTEGER NOT NULL,
  name VARCHAR(50) NOT NULL,
  description VARCHAR(10000) NOT NULL
);
-- ---
-- Table ‘nmap_network’
-- 
-- ---
DROP TABLE IF EXISTS nmap_network;
CREATE TABLE nmap_network (
  id serial primary key,
  name VARCHAR(500) NOT NULL
);
-- ---
-- Table ‘cellline’
-- 
-- ---
DROP TABLE IF EXISTS cellline;
CREATE TABLE cellline (
  id serial primary key,
  name VARCHAR(500) NOT NULL,
  description VARCHAR(10000)
);
-- ---
-- Table ‘nmap_node’
-- 
-- ---
DROP TABLE IF EXISTS nmap_node;
CREATE TABLE nmap_node (
  id serial primary key,
  type_id INTEGER NOT NULL,
  name VARCHAR(500) NOT NULL
);
-- ---
-- Table ‘nmap_edge’
-- 
-- ---
DROP TABLE IF EXISTS nmap_edge;
CREATE TABLE nmap_edge (
  id serial primary key,
  start INTEGER NOT NULL,
  "end" INTEGER NOT NULL,
  model_id INTEGER NOT NULL,
  edge_activity BOOLEAN NOT NULL
);
-- ---
-- Table ‘node_type’
-- 
-- ---
DROP TABLE IF EXISTS node_type;
CREATE TABLE node_type (
  id serial primary key,
  name VARCHAR(500) NOT NULL
);
-- ---
-- Table ‘node_activity’
-- 
-- ---
DROP TABLE IF EXISTS node_activity;
CREATE TABLE node_activity (
  id serial primary key,
  network_id INTEGER NOT NULL,
  node_id INTEGER NOT NULL,
  activity_id INTEGER NOT NULL
);
-- ---
-- Table ‘node_perturbation’
-- 
-- ---
DROP TABLE IF EXISTS node_perturbation;
CREATE TABLE node_perturbation (
  id serial primary key,
  node_id INTEGER NOT NULL,
  activity_id INTEGER NOT NULL,
  perturbation_id INTEGER NOT NULL
);
-- ---
-- Table ‘node_activity_type’
-- 
-- ---
DROP TABLE IF EXISTS node_activity_type;
CREATE TABLE node_activity_type (
  id serial primary key,
  name VARCHAR(100) NOT NULL
);
-- ---
-- Table ‘perturbation’
-- 
-- ---
DROP TABLE IF EXISTS perturbation;
CREATE TABLE perturbation (
  id serial primary key,
  network_id INTEGER NOT NULL,
  name VARCHAR(500) NOT NULL,
  description VARCHAR(10000) NOT NULL
);
-- ---
-- Table ‘attractor’
-- 
-- ---
DROP TABLE IF EXISTS attractor;
CREATE TABLE attractor (
  id serial primary key,
  perturbation_id INTEGER NOT NULL,
  ratio FLOAT NOT NULL
);
-- ---
-- Table ‘nmap_phenotype’
-- 
-- ---
DROP TABLE IF EXISTS nmap_phenotype;
CREATE TABLE nmap_phenotype (
  id serial primary key,
  name VARCHAR(500) NOT NULL
);
-- ---
-- Table ‘attractor_phenotype’
-- 
-- ---
DROP TABLE IF EXISTS attractor_phenotype;
CREATE TABLE attractor_phenotype (
  id serial primary key,
  attractor_id INTEGER NOT NULL,
  phenotype_id INTEGER NOT NULL,
  ratio FLOAT NOT NULL
);
-- ---
-- Table ‘network_cellline’
-- 
-- ---
DROP TABLE IF EXISTS network_cellline;
CREATE TABLE network_cellline (
  id serial primary key,
  network_id INTEGER NOT NULL,
  cellline_id INTEGER NOT NULL
);

-- ---
-- Table ‘taget_node’
-- 
-- ---
-- DROP TABLE IF EXISTS target_node;
-- CREATE TABLE target_node( 
--   id bigserial primary key, 
--   protein_id BIGINT NOT NULL, 
--   node_id BIGINT NOT NULL
-- );

--

-- ---
-- Table 'tag_table'
--
-- ---
DROP TABLE IF EXISTS tag_table;
CREATE TABLE tag_table (
  id BIGINT NOT NULL PRIMARY KEY
);

--

-- ---
-- Table 'physchem'
--
-- ---
--DROP TABLE IF EXISTS physchem;
--CREATE TABLE physchem (
-- id BIGINT primary key,
-- mw FLOAT NOT NULL,
-- logp FLOAT NULL DEFAULT NULL,
-- tpsa FLOAT NULL DEFAULT NULL,
-- mw_class INTEGER NOT NULL,
-- logp_class INTEGER NOT NULL,
-- solubility_class INTEGER NOT NULL
--);

-- ---
-- Table 'lipinski'
--
-- ---
--DROP TABLE IF EXISTS lipinski;
--CREATE TABLE lipinski (
-- id BIGINT primary key,
-- h_donors INTEGER NULL DEFAULT NULL,
-- h_acceptors INTEGER NULL DEFAULT NULL,
-- rot_bonds INTEGER NULL DEFAULT NULL,
-- rings INTEGER NULL DEFAULT NULL,
-- h_donors_class INTEGER NOT NULL,
-- h_acceptors_class INTEGER NOT NULL,
-- rot_bonds_class INTEGER NOT NULL,
-- rings_class INTEGER NOT NULL
--);

-- ---
-- Table 'admet'
--
-- ---
--DROP TABLE IF EXISTS admet;
--CREATE TABLE admet (
-- id BIGINT primary key,
-- caco2_permeability_class INTEGER NULL DEFAULT NULL,
-- hia_class INTEGER NULL DEFAULT NULL,
-- metabolic_stability_class INTEGER NOT NULL,
-- ames_test_class INTEGER NOT NULL,
-- herg_inhibition_class INTEGER NOT NULL
--);

-- ---
-- Table 'mw_class_type'
--
-- ---
DROP TABLE IF EXISTS mw_class_type;
CREATE TABLE mw_class_type (
 id INTEGER NOT NULL,
 name VARCHAR(20) NOT NULL,
 PRIMARY KEY (id)
);

-- ---
-- Table 'logp_class_type'
--
-- ---
DROP TABLE IF EXISTS logp_class_type;
CREATE TABLE logp_class_type (
 id INTEGER primary key,
 name VARCHAR(20) NOT NULL
);

-- ---
-- Table 'solubility_class_type'
--
-- ---
DROP TABLE IF EXISTS solubility_class_type;
CREATE TABLE solubility_class_type (
 id INTEGER primary key,
 name VARCHAR(20) NOT NULL
);

-- ---
-- Table 'lipinski_class_type'
--
-- ---
-- DROP TABLE IF EXISTS lipinski_class_type;
-- CREATE TABLE lipinski_class_type (
--  id INTEGER primary key,
--  name VARCHAR(20) NOT NULL
-- );

-- ---
-- Table 'permeability_class_type'
--
-- ---
DROP TABLE IF EXISTS permeability_class_type;
CREATE TABLE permeability_class_type (
 id INTEGER primary key,
 name VARCHAR(20) NOT NULL
);

-- ---
-- Table 'hia_class_type'
--
-- ---
DROP TABLE IF EXISTS hia_class_type;
CREATE TABLE hia_class_type (
 id INTEGER primary key,
 name VARCHAR(20) NOT NULL
);

-- ---
-- Table 'stability_class_type'
--
-- ---
DROP TABLE IF EXISTS stability_class_type;
CREATE TABLE stability_class_type (
 id INTEGER primary key,
 name VARCHAR(20) NOT NULL
);

-- ---
-- Table 'ames_class_type'
--
-- ---
DROP TABLE IF EXISTS ames_class_type;
CREATE TABLE ames_class_type (
 id INTEGER primary key,
 name VARCHAR(20) NOT NULL
);

-- ---
-- Table 'herg_class_type'
--
-- ---
DROP TABLE IF EXISTS herg_class_type;
CREATE TABLE herg_class_type (
 id INTEGER primary key,
 name VARCHAR(20) NOT NULL
);

-- ---
-- Table 'feature_fingerprint'
--
-- ---

DROP TABLE IF EXISTS feature_fingerprint;
CREATE TABLE feature_fingerprint (
   id bigserial primary key,
  conformer_id BIGINT NOT NULL,
  type CHAR(3) NOT NULL,
  fp_triad BIGINT NOT NULL
);

-- ---
-- Table 'atom3d'
--
-- ---

DROP TABLE IF EXISTS atom3d;
CREATE TABLE atom3d (
   id bigserial primary key,
  conformer_id BIGINT NOT NULL,
  element_id INTEGER NOT NULL,
  x FLOAT NOT NULL,
  y FLOAT NOT NULL,
  z FLOAT NOT NULL
);

-- ---
-- Table 'element'
--
-- ---

DROP TABLE IF EXISTS element;
CREATE TABLE element (
   id serial primary key,
  name VARCHAR(20) NOT NULL,
  symbol VARCHAR(2)  NOT NULL,
  number INTEGER NOT NULL
);

-- ---
-- Table 'env_all'
--
-- ---

DROP TABLE IF EXISTS env_all;
CREATE TABLE env_all (
  id BIGINT primary key,
  freq00 FLOAT NOT NULL DEFAULT 0,
  freq01 FLOAT NOT NULL DEFAULT 0,
  freq02 FLOAT NOT NULL DEFAULT 0,
  freq03 FLOAT NOT NULL DEFAULT 0,
  freq04 FLOAT NOT NULL DEFAULT 0,
  freq05 FLOAT NOT NULL DEFAULT 0,
  freq06 FLOAT NOT NULL DEFAULT 0
);

-- ---
-- Table 'env_c'
--
-- ---

DROP TABLE IF EXISTS env_c;
CREATE TABLE env_c (
  id BIGINT primary key,
  freq00 FLOAT NOT NULL DEFAULT 0,
  freq01 FLOAT NOT NULL DEFAULT 0,
  freq02 FLOAT NOT NULL DEFAULT 0,
  freq03 FLOAT NOT NULL DEFAULT 0,
  freq04 FLOAT NOT NULL DEFAULT 0,
  freq05 FLOAT NOT NULL DEFAULT 0,
  freq06 FLOAT NOT NULL DEFAULT 0
);

-- ---
-- Table 'env_n'
--
-- ---

DROP TABLE IF EXISTS env_n;
CREATE TABLE env_n (
  id BIGINT primary key,
  freq00 FLOAT NOT NULL DEFAULT 0,
  freq01 FLOAT NOT NULL DEFAULT 0,
  freq02 FLOAT NOT NULL DEFAULT 0,
  freq03 FLOAT NOT NULL DEFAULT 0,
  freq04 FLOAT NOT NULL DEFAULT 0,
  freq05 FLOAT NOT NULL DEFAULT 0,
  freq06 FLOAT NOT NULL DEFAULT 0
);

-- ---
-- Table 'env_o'
--
-- ---

DROP TABLE IF EXISTS env_o;
CREATE TABLE env_o (
  id BIGINT primary key,
  freq00 FLOAT NOT NULL DEFAULT 0,
  freq01 FLOAT NOT NULL DEFAULT 0,
  freq02 FLOAT NOT NULL DEFAULT 0,
  freq03 FLOAT NOT NULL DEFAULT 0,
  freq04 FLOAT NOT NULL DEFAULT 0,
  freq05 FLOAT NOT NULL DEFAULT 0,
  freq06 FLOAT NOT NULL DEFAULT 0
);

-- ---
-- Table 'env_p'
--
-- ---

DROP TABLE IF EXISTS env_p;
CREATE TABLE env_p (
  id BIGINT primary key,
  freq00 FLOAT NOT NULL DEFAULT 0,
  freq01 FLOAT NOT NULL DEFAULT 0,
  freq02 FLOAT NOT NULL DEFAULT 0,
  freq03 FLOAT NOT NULL DEFAULT 0,
  freq04 FLOAT NOT NULL DEFAULT 0,
  freq05 FLOAT NOT NULL DEFAULT 0,
  freq06 FLOAT NOT NULL DEFAULT 0
);

-- ---
-- Table 'env_s'
--
-- ---

DROP TABLE IF EXISTS env_s;
CREATE TABLE env_s (
  id BIGINT primary key,
  freq00 FLOAT NOT NULL DEFAULT 0,
  freq01 FLOAT NOT NULL DEFAULT 0,
  freq02 FLOAT NOT NULL DEFAULT 0,
  freq03 FLOAT NOT NULL DEFAULT 0,
  freq04 FLOAT NOT NULL DEFAULT 0,
  freq05 FLOAT NOT NULL DEFAULT 0,
  freq06 FLOAT NOT NULL DEFAULT 0
);

-- ---
-- Table 'env_halogen'
--
-- ---

DROP TABLE IF EXISTS env_halogen;
CREATE TABLE env_halogen (
  id BIGINT primary key,
  freq00 FLOAT NOT NULL DEFAULT 0,
  freq01 FLOAT NOT NULL DEFAULT 0,
  freq02 FLOAT NOT NULL DEFAULT 0,
  freq03 FLOAT NOT NULL DEFAULT 0,
  freq04 FLOAT NOT NULL DEFAULT 0,
  freq05 FLOAT NOT NULL DEFAULT 0,
  freq06 FLOAT NOT NULL DEFAULT 0
);

-- ---
-- Table 'env_fingerprint'
--
-- ---

DROP TABLE IF EXISTS env_fingerprint;
CREATE TABLE env_fingerprint (
  id BIGINT primary key,
  fp_all_01 BIGINT NOT NULL DEFAULT 0,
  fp_all_02 BIGINT NOT NULL DEFAULT 0,
  fp_all_03 BIGINT NOT NULL DEFAULT 0,
  fp_c_01 BIGINT NOT NULL DEFAULT 0,
  fp_c_02 BIGINT NOT NULL DEFAULT 0,
  fp_c_03 BIGINT NOT NULL DEFAULT 0,
  fp_n_01 BIGINT NOT NULL DEFAULT 0,
  fp_n_02 BIGINT NOT NULL DEFAULT 0,
  fp_n_03 BIGINT NOT NULL DEFAULT 0,
  fp_o_01 BIGINT NOT NULL DEFAULT 0,
  fp_o_02 BIGINT NOT NULL DEFAULT 0,
  fp_o_03 BIGINT NOT NULL DEFAULT 0,
  fp_p_01 BIGINT NOT NULL DEFAULT 0,
  fp_p_02 BIGINT NOT NULL DEFAULT 0,
  fp_p_03 BIGINT NOT NULL DEFAULT 0,
  fp_s_01 BIGINT NOT NULL DEFAULT 0,
  fp_s_02 BIGINT NOT NULL DEFAULT 0,
  fp_s_03 BIGINT NOT NULL DEFAULT 0,
  fp_hal_01 BIGINT NOT NULL DEFAULT 0,
  fp_hal_02 BIGINT NOT NULL DEFAULT 0,
  fp_hal_03 BIGINT NOT NULL DEFAULT 0
);

-- ---
-- Table 'physchem'
--
-- ---

DROP TABLE IF EXISTS physchem;
CREATE TABLE physchem (
  id BIGINT primary key,
  mw FLOAT NOT NULL,
  logp FLOAT NULL DEFAULT NULL,
  tpsa FLOAT NULL DEFAULT NULL,
  h_donors INTEGER NULL DEFAULT NULL,
  h_acceptors INTEGER NULL DEFAULT NULL,
  rot_bonds INTEGER NULL DEFAULT NULL,
  rings INTEGER NULL DEFAULT NULL
);

-- ---
-- Table 'absorption'
--
-- ---

DROP TABLE IF EXISTS absorption;
CREATE TABLE absorption (
  id BIGINT primary key,
  a FLOAT NULL DEFAULT NULL,
  b FLOAT NULL DEFAULT NULL,
  bo FLOAT NULL DEFAULT NULL,
  l FLOAT NULL DEFAULT NULL,
  s FLOAT NULL DEFAULT NULL,
  e FLOAT NULL DEFAULT NULL,
  mcgowanvolume FLOAT NULL DEFAULT NULL,
  abs_oral_dose_02mg FLOAT NULL DEFAULT NULL,
  abs_oral_dose_1mg FLOAT NULL DEFAULT NULL,
  abs_oral_dose_5mg FLOAT NULL DEFAULT NULL,
  abs_oral_dose_20mg FLOAT NULL DEFAULT NULL,
  abs_oral_dose_50mg FLOAT NULL DEFAULT NULL,
  abs_oral_dose_100mg FLOAT NULL DEFAULT NULL,
  abs_oral_dose_500mg FLOAT NULL DEFAULT NULL,
  pe FLOAT NULL DEFAULT NULL,
  ka FLOAT NULL DEFAULT NULL,
  maxhia FLOAT NULL DEFAULT NULL,
  pecaco2_ph65_0rpm FLOAT NULL DEFAULT NULL,
  pecaco2_ph65_500rpm FLOAT NULL DEFAULT NULL,
  pecaco2_ph74_0rpm FLOAT NULL DEFAULT NULL,
  pecaco2_ph74_500rpm FLOAT NULL DEFAULT NULL,
  pecaco2_ph80_0rpm FLOAT NULL DEFAULT NULL,
  pecaco2_ph80_500rpm FLOAT NULL DEFAULT NULL,
  mrdd FLOAT NULL DEFAULT NULL,
  mrddclass FLOAT NULL DEFAULT NULL
);

-- ---
-- Table 'distribution'
--
-- ---

DROP TABLE IF EXISTS distribution;
CREATE TABLE distribution (
  id BIGINT primary key,
  logps FLOAT NULL DEFAULT NULL,
  logbb FLOAT NULL DEFAULT NULL,
  fu_brain FLOAT NULL DEFAULT NULL,
  bbb_eq_rate FLOAT NULL DEFAULT NULL,
  _ppb FLOAT NULL DEFAULT NULL,
  _ppb_ri FLOAT NULL DEFAULT NULL,
  logka_hsa FLOAT NULL DEFAULT NULL,
  logka_hsa_ri FLOAT NULL DEFAULT NULL,
  vss FLOAT NULL DEFAULT NULL,
  actr_asbt_substrate FLOAT NULL DEFAULT NULL,
  actr_pept1_substrate FLOAT NULL DEFAULT NULL,
  actr_mct1_substrate FLOAT NULL DEFAULT NULL,
  actr_amino_acid_transport FLOAT NULL DEFAULT NULL,
  actr_other_transport FLOAT NULL DEFAULT NULL
);

-- ---
-- Table 'metabolism'
--
-- ---

DROP TABLE IF EXISTS metabolism;
CREATE TABLE metabolism (
  id BIGINT primary key,
  _pgp_s FLOAT NULL DEFAULT NULL,
  _pgp_s_ri FLOAT NULL DEFAULT NULL,
  _pgp_s_good FLOAT NULL DEFAULT NULL,
  _pgp_s_good_ri FLOAT NULL DEFAULT NULL,
  _pgp_i FLOAT NULL DEFAULT NULL,
  _pgp_i_ri FLOAT NULL DEFAULT NULL,
  _pgp_i_good FLOAT NULL DEFAULT NULL,
  _pgp_i_good_ri FLOAT NULL DEFAULT NULL,
  _cyp3a4_i_10 FLOAT NULL DEFAULT NULL,
  _cyp3a4_i_10_ri FLOAT NULL DEFAULT NULL,
  _cyp2d6_i_10 FLOAT NULL DEFAULT NULL,
  _cyp2d6_i_10_ri FLOAT NULL DEFAULT NULL,
  _cyp2c9_i_10 FLOAT NULL DEFAULT NULL,
  _cyp2c9_i_10_ri FLOAT NULL DEFAULT NULL,
  _cyp2c19_i_10 FLOAT NULL DEFAULT NULL,
  _cyp2c19_i_10_ri FLOAT NULL DEFAULT NULL,
  _cyp1a2_i_10 FLOAT NULL DEFAULT NULL,
  _cyp1a2_i_10_ri FLOAT NULL DEFAULT NULL,
  _cyp3a4_i_50 FLOAT NULL DEFAULT NULL,
  _cyp3a4_i_50_ri FLOAT NULL DEFAULT NULL,
  _cyp2d6_i_50 FLOAT NULL DEFAULT NULL,
  _cyp2d6_i_50_ri FLOAT NULL DEFAULT NULL,
  _cyp2c9_i_50 FLOAT NULL DEFAULT NULL,
  _cyp2c9_i_50_ri FLOAT NULL DEFAULT NULL,
  _cyp2c19_i_50 FLOAT NULL DEFAULT NULL,
  _cyp2c19_i_50_ri FLOAT NULL DEFAULT NULL,
  _cyp1a2_i_50 FLOAT NULL DEFAULT NULL,
  _cyp1a2_i_50_ri FLOAT NULL DEFAULT NULL,
  _cyp3a4_s FLOAT NULL DEFAULT NULL,
  _cyp3a4_s_ri FLOAT NULL DEFAULT NULL,
  _cyp2d6_s FLOAT NULL DEFAULT NULL,
  _cyp2d6_s_ri FLOAT NULL DEFAULT NULL,
  _cyp2c9_s FLOAT NULL DEFAULT NULL,
  _cyp2c9_s_ri FLOAT NULL DEFAULT NULL,
  _cyp2c19_s FLOAT NULL DEFAULT NULL,
  _cyp2c19_s_ri FLOAT NULL DEFAULT NULL,
  _cyp1a2_s FLOAT NULL DEFAULT NULL,
  _cyp1a2_s_ri FLOAT NULL DEFAULT NULL
);

-- ---
-- Table 'toxicity'
--
-- ---

DROP TABLE IF EXISTS toxicity;
CREATE TABLE toxicity (
  id BIGINT primary key,
  _ld50_mouse_ip FLOAT NULL DEFAULT NULL,
  _ld50_mouse_ip_ri FLOAT NULL DEFAULT NULL,
  _ld50_mouse_or FLOAT NULL DEFAULT NULL,
  _ld50_mouse_or_ri FLOAT NULL DEFAULT NULL,
  _ld50_mouse_iv FLOAT NULL DEFAULT NULL,
  _ld50_mouse_iv_ri FLOAT NULL DEFAULT NULL,
  _ld50_mouse_sc FLOAT NULL DEFAULT NULL,
  _ld50_mouse_sc_ri FLOAT NULL DEFAULT NULL,
  _ld50_rat_ip FLOAT NULL DEFAULT NULL,
  _ld50_rat_ip_ri FLOAT NULL DEFAULT NULL,
  _ld50_rat_or FLOAT NULL DEFAULT NULL,
  _ld50_rat_or_ri FLOAT NULL DEFAULT NULL,
  _lc50_ppromelas FLOAT NULL DEFAULT NULL,
  _lc50_ppromelas_ri FLOAT NULL DEFAULT NULL,
  _lc50_dmagna FLOAT NULL DEFAULT NULL,
  _lc50_dmagna_ri FLOAT NULL DEFAULT NULL,
  _igc50_tpyriformis FLOAT NULL DEFAULT NULL,
  _igc50_tpyriformis_ri FLOAT NULL DEFAULT NULL,
  _ames FLOAT NULL DEFAULT NULL,
  _ames_ri FLOAT NULL DEFAULT NULL,
  _logrba_3 FLOAT NULL DEFAULT NULL,
  _logrba_3_ri FLOAT NULL DEFAULT NULL,
  _logrba_0 FLOAT NULL DEFAULT NULL,
  _logrba_0_ri FLOAT NULL DEFAULT NULL,
  probability_of_he_on_blood FLOAT NULL DEFAULT NULL,
  probability_of_he_on_cardiovascular_system FLOAT NULL DEFAULT NULL,
  probability_of_he_on_gastrointestinal_system FLOAT NULL DEFAULT NULL,
  probability_of_he_on_kidney FLOAT NULL DEFAULT NULL,
  probability_of_he_on_liver FLOAT NULL DEFAULT NULL,
  probability_of_he_on_lungs FLOAT NULL DEFAULT NULL,
  _herg FLOAT NULL DEFAULT NULL,
  _herg_ri FLOAT NULL DEFAULT NULL,
  eye_irritation FLOAT NULL DEFAULT NULL,
  skin_irritation FLOAT NULL DEFAULT NULL
);

-- ---
-- Table 'leadlike_category'
--
-- ---

DROP TABLE IF EXISTS leadlike_category;
CREATE TABLE leadlike_category (
  id BIGINT primary key,
  mw_class INTEGER NOT NULL,
  logp_class INTEGER NOT NULL,
  solubility_class INTEGER NOT NULL,
  caco2_permeability_class INTEGER NULL DEFAULT NULL,
  hia_class INTEGER NULL DEFAULT NULL,
  metabolic_stability_class INTEGER NOT NULL,
  ames_test_class INTEGER NOT NULL,
  herg_inhibition_class INTEGER NOT NULL
);



-- ---
-- Table 'gene'
--
-- ---
DROP TABLE IF EXISTS gene;
CREATE TABLE gene (
 id serial primary key,
 ensembl_id VARCHAR(30) NOT NULL,
 name VARCHAR(30) NULL DEFAULT NULL,
 description VARCHAR(200) NULL DEFAULT NULL,
 loc_type VARCHAR(30) NOT NULL,
 loc_seq VARCHAR(50) NOT NULL,
 start INTEGER NOT NULL,
 "end" INTEGER NOT NULL
);


-- ---
-- Table 'hgnc'
--
-- ---
DROP TABLE IF EXISTS hgnc;
CREATE TABLE hgnc (
 id serial primary key,
 hgnc_id VARCHAR(20) NOT NULL,
 gene_id INTEGER NOT NULL
);


-- ---
-- Table 'gene_node'
--
-- ---
DROP TABLE IF EXISTS gene_node;
CREATE TABLE gene_node (
 id serial primary key,
 gene_id INTEGER NOT NULL,
 node_id INTEGER NOT NULL
);


-- ---
-- Table 'target_molecule'
--
-- ---
DROP TABLE IF EXISTS target_molecule;
CREATE TABLE target_molecule (
  id serial primary key,
  protein_id BIGINT NOT NULL,
  molecule_id BIGINT NOT NULL,
  activity FLOAT NOT NULL
);

-- ---
-- Table 'protein_chembl'
--
-- ---
DROP TABLE IF EXISTS protein_chembl;
CREATE TABLE protein_chembl (
  id serial primary key,
  protein_id BIGINT NOT NULL,
  chembl_id VARCHAR(20) NOT NULL
);

-- ---
-- Ntropy utils made in onmakers
-- ---

-- ---
-- Table 'tc_access_log'
--
-- ---
DROP TABLE IF EXISTS tc_access_log;
CREATE TABLE tc_access_log (
	log_seq serial   NOT null PRIMARY KEY,
	rqst_url varchar(255) NULL,
	query_string varchar(1000) NULL,
	referer varchar(1000) NULL,
	cont_ip varchar(30) NULL,
	locale varchar(20) NULL,
	user_agent varchar(255) NULL,
	device_divn varchar(20) NULL,
	os_nm varchar(20) NULL,
	browser_nm varchar(20) NULL,
	browser_vrsn varchar(20) NULL,
	session_id varchar(100) NULL,
	resolution varchar(20) NULL,
	creat_dt timestamptz NULL,
	menu_key varchar(255) NULL,
	user_id varchar(30) null
);

-- ---
-- Table 'tc_access_stats_brwsr'
--
-- ---
DROP TABLE IF EXISTS tc_access_stats_brwsr;
CREATE TABLE tc_access_stats_brwsr (
	id serial NOT null PRIMARY KEY,
	stts_dt varchar(10) NULL,
	browser_nm varchar(20) NULL,
	browser_vrsn varchar(20) NULL,
	device_divn varchar(20) NULL,
	vsit_cnt int4 NOT NULL,
	page_read_cnt int4 NOT NULL
);

-- ---
-- Table 'tc_access_stats_de'
--
-- ---
DROP TABLE IF EXISTS tc_access_stats_de;
CREATE TABLE tc_access_stats_de (
	id serial NOT null PRIMARY KEY,
	stts_dt varchar(10) NULL,
	vsit_cnt int4 NOT NULL,
	user_vsit_cnt int4 NOT NULL,
	guest_vsit_cnt int4 NOT NULL,
	page_read_cnt int4 NOT NULL,
	user_page_read_cnt int4 NOT NULL,
	guest_page_read_cnt int4 NOT NULL
);

-- ---
-- Table 'tc_access_stats_menu'
--
-- ---
DROP TABLE IF EXISTS tc_access_stats_menu;
CREATE TABLE tc_access_stats_menu (
	id serial NOT null PRIMARY KEY,
	stts_dt varchar(10) NULL,
	menu_key varchar(36) NULL,
	vsit_cnt int4 NOT NULL,
	user_vsit_cnt int4 NOT NULL,
	guest_vsit_cnt int4 NOT NULL,
	page_read_cnt int4 NOT NULL,
	user_page_read_cnt int4 NOT NULL,
	guest_page_read_cnt int4 NOT NULL
);

-- ---
-- Table 'tc_access_stats_opersysm'
--
-- ---
DROP TABLE IF EXISTS tc_access_stats_opersysm;
CREATE TABLE tc_access_stats_opersysm (
	id serial NOT null PRIMARY KEY,
	stts_dt varchar(10) NULL,
	os_nm varchar(20) NULL,
	device_divn varchar(20) NULL,
	vsit_cnt int4 NOT NULL,
	page_read_cnt int4 NOT NULL
);

-- ---
-- Table 'tc_access_stats_rsoltn'
--
-- ---
DROP TABLE IF EXISTS tc_access_stats_rsoltn;
CREATE TABLE tc_access_stats_rsoltn (
	id serial NOT null PRIMARY KEY,
	stts_dt varchar(10) NULL,
	resolution varchar(20) NULL,
	device_divn varchar(20) NULL,
	vsit_cnt int4 NOT NULL,
	page_read_cnt int4 NOT NULL
);

-- ---
-- Table 'tc_banner'
--
-- ---
DROP TABLE IF EXISTS tc_banner;
CREATE TABLE tc_banner (
	banner_key varchar(36) NULL,
	banner_title varchar(2000) NULL,
	ntce_bgnde date NULL,
	ntce_endde date NULL,
	link_adres varchar(300) NULL,
	attach_nm varchar(1000) NULL,
	attach_size int4 NOT NULL,
	attach_save_nm varchar(50) NULL,
	attach_save_path varchar(200) NULL,
	act_yn varchar(1) NULL,
	use_yn varchar(1) NULL,
	sort_sn int4 NOT NULL,
	creat_user varchar(100) NULL,
	creat_dt timestamptz NULL,
	updt_user varchar(100) NULL,
	updt_dt timestamptz NULL
);

-- ---
-- Table 'tc_bbs'
--
-- ---
DROP TABLE IF EXISTS tc_bbs;
CREATE TABLE tc_bbs (
	board_key varchar(36) NULL primary key,
	board_nm varchar(500) NULL,
	ctgry_grp varchar(20) NULL,
	board_type_cd varchar(20) NULL,
	use_yn varchar(1) NULL,
	notification_fn_yn varchar(1) NULL,
	comment_fn_yn varchar(1) NULL,
	private_fn_yn varchar(1) NULL,
	attach_fn_yn varchar(1) NULL,
	attach_permit_cnt int4 NOT NULL,
	new_notice_term_days int4 NOT NULL,
	attach_permit_size int4 NOT NULL,
	page_post_cnt int4 NOT NULL,
	creat_user varchar(100) NULL,
	creat_dt timestamptz NULL,
	updt_user varchar(100) NULL,
	updt_dt timestamptz NULL
);

-- ---
-- Table 'tc_bbsctt'
--
-- ---
DROP TABLE IF EXISTS tc_bbsctt;
CREATE TABLE tc_bbsctt (
	post_key varchar(36) NULL primary key,
	title varchar(2000) NULL,
	hit_cnt int4 NOT NULL,
	post_close_yn varchar(1) NULL,
	notification_yn varchar(1) NULL,
	post_passwd varchar(150) NULL,
	use_yn varchar(1) NULL,
	reg_dt date NULL,
	reg_user_nm varchar(100) NULL,
	creat_user varchar(100) NULL,
	creat_dt timestamptz NULL,
	updt_user varchar(100) NULL,
	updt_dt timestamptz NULL,
	board_key varchar(36) NULL
);

-- ---
-- Table 'tc_bbsctt_atchmnfl'
--
-- ---
DROP TABLE IF EXISTS tc_bbsctt_atchmnfl;
CREATE TABLE tc_bbsctt_atchmnfl (
	attach_key varchar(36) NOT NULL,
	attach_nm varchar(1000) NULL,
	attach_size int4 NOT NULL,
	attach_save_nm varchar(100) NULL,
	attach_save_path varchar(200) NULL,
	content_type_cd varchar(20) NULL,
	use_yn varchar(1) NULL,
	creat_user varchar(100) NULL,
	creat_dt timestamptz NULL,
	updt_user varchar(100) NULL,
	updt_dt timestamptz NULL,
	post_key varchar(36) NULL,
	CONSTRAINT tc_bbsctt_atchmnfl_pkey PRIMARY KEY (attach_key)
);

-- ---
-- Table 'tc_bbsctt_cn'
--
-- ---
DROP TABLE IF EXISTS tc_bbsctt_cn;
CREATE TABLE tc_bbsctt_cn (
	id serial NOT null PRIMARY KEY,
	content_type_cd varchar(20) NULL,
	"content" text NULL,
	reg_reply_user_nm varchar(100) NULL,
	creat_user varchar(100) NULL,
	creat_dt timestamptz NULL,
	updt_user varchar(100) NULL,
	updt_dt timestamptz NULL,
	post_key varchar(36) NULL
);

-- ---
-- Table 'tc_bbsctt_cn_log'
--
-- ---
DROP TABLE IF EXISTS tc_bbsctt_cn_log;
CREATE TABLE tc_bbsctt_cn_log (
	id serial NOT null PRIMARY KEY,
	post_key varchar(36) NULL,
	content_type_cd varchar(20) NULL,
	"content" text NULL,
	reg_reply_user_nm varchar(100) NULL,
	creat_user varchar(100) NULL,
	creat_dt timestamptz NULL,
	updt_user varchar(100) NULL,
	updt_dt timestamptz NULL
);

-- ---
-- Table 'tc_bbsctt_comment'
--
-- ---
DROP TABLE IF EXISTS tc_bbsctt_comment;
CREATE TABLE tc_bbsctt_comment (
	comment_key varchar(36) NULL,
	comment_password varchar(1000) NULL,
	"comment" text NULL,
	use_yn varchar(1) NULL,
	reg_comment_user_nm varchar(100) NULL,
	creat_user varchar(100) NULL,
	creat_dt timestamptz NULL,
	updt_user varchar(100) NULL,
	updt_dt timestamptz NULL,
	post_key varchar(36) NULL
);

-- ---
-- Table 'tc_cd'
--
-- ---
DROP TABLE IF EXISTS tc_cd;
CREATE TABLE tc_cd (
	code varchar(20) NULL,
	code_grp varchar(20) NULL,
	code_nm varchar(500) NULL,
	use_yn varchar(1) NULL,
	creat_user varchar(100) NULL,
	creat_dt timestamptz NULL,
	updt_user varchar(100) NULL,
	updt_dt timestamptz NULL
);

-- ---
-- Table 'tc_menu'
--
-- ---
DROP TABLE IF EXISTS tc_menu;
CREATE TABLE tc_menu (
	menu_key varchar(36) NULL primary key,
	up_menu_key varchar(36) NULL,
	menu_nm varchar(500) NULL,
	menu_type_cd varchar(20) NULL,
	menu_sort int4 NOT NULL,
	link_url varchar(2000) NULL,
	link_new_window_yn varchar(1) NULL,
	post_key varchar(36) NULL,
	layout varchar(50) NULL,
	board_key varchar(36) NULL,
	menu_desc varchar(2000) NULL,
	img_path varchar(100) NULL,
	use_yn varchar(1) NULL,
	creat_user varchar(100) NULL,
	creat_dt timestamptz NULL,
	updt_dt timestamptz NULL,
	updt_user varchar(100) NULL,
	access_mng varchar(1) NULL,
	role_id varchar(30) NULL
);

-- ---
-- Table 'tc_menu_mngr'
--
-- ---
DROP TABLE IF EXISTS tc_menu_mngr;
CREATE TABLE tc_menu_mngr (
	id serial NOT null PRIMARY KEY,
	edit_auth varchar(1) NULL,
	creat_user varchar(100) NULL,
	creat_dt timestamptz NULL,
	updt_user varchar(100) NULL,
	updt_dt timestamptz NULL,
	menu_key varchar(36) NULL,
	member_id int4 NULL
);

-- ---
-- Table 'tc_mngr'
--
-- ---
DROP TABLE IF EXISTS tc_mngr;
CREATE TABLE tc_mngr (
	admin_id varchar(20) NULL,
	admin_nm varchar(50) NULL,
	admin_password varchar(150) NULL,
	use_yn varchar(1) NULL,
	creat_user varchar(100) NULL,
	creat_dt timestamptz NULL,
	updt_user varchar(100) NULL,
	updt_dt timestamptz NULL
);

-- ---
-- Table 'tc_popup'
--
-- ---
DROP TABLE IF EXISTS tc_popup;
CREATE TABLE tc_popup (
	popup_key varchar(36) NULL,
	popup_title varchar(2000) NULL,
	ntce_bgnde date NULL,
	ntce_endde date NULL,
	link_adres varchar(300) NULL,
	width_size int4 NOT NULL,
	height_size int4 NOT NULL,
	top_position int4 NOT NULL,
	left_position int4 NOT NULL,
	"content" text NULL,
	act_yn varchar(1) NULL,
	use_yn varchar(1) NULL,
	sort_sn int4 NOT NULL,
	creat_user varchar(100) NULL,
	creat_dt timestamptz NULL,
	updt_user varchar(100) NULL,
	updt_dt timestamptz NULL
);

-- ---
-- Table 'tc_vis'
--
-- ---
DROP TABLE IF EXISTS tc_vis;
CREATE TABLE tc_vis (
	visual_key varchar(36) NULL,
	visual_title varchar(2000) NULL,
	visual_sub_title varchar(2000) NULL,
	ntce_bgnde date NULL,
	ntce_endde date NULL,
	link_adres varchar(300) NULL,
	attach_nm varchar(1000) NULL,
	attach_size int4 NOT NULL,
	attach_save_nm varchar(50) NULL,
	attach_save_path varchar(200) NULL,
	act_yn varchar(1) NULL,
	use_yn varchar(1) NULL,
	sort_sn int4 NOT NULL,
	creat_user varchar(100) NULL,
	creat_dt timestamptz NULL,
	updt_user varchar(100) NULL,
	updt_dt timestamptz NULL
);

-- ---
-- Table 'member_role'
--
-- ---
DROP TABLE IF EXISTS member_role;
CREATE TABLE member_role (
	role_id varchar(30) NULL primary key,
	role_name varchar(100) NULL,
	role_level int4 NOT NULL,
	creat_user varchar(100) NULL,
	creat_dt timestamptz NULL,
	updt_user varchar(100) NULL,
	updt_dt timestamptz NULL,
	use_yn bpchar(1) NULL,
	sort_sn int4 NULL
);

-- ---
-- Table 'project_member_role'
--
-- ---
DROP TABLE IF EXISTS project_member_role;
CREATE TABLE project_member_role (
	project_role_id varchar(30) NOT NULL,
	project_role_name varchar(100) NULL,
	project_role_level int4 NOT NULL,
	use_yn bpchar(1) NULL,
	creat_user varchar(30) NULL,
	creat_dt timestamp NULL,
	updt_user varchar(30) NULL,
	updt_dt timestamp NULL,
	CONSTRAINT pk_project_member_role PRIMARY KEY (project_role_id)
);

-- ---
-- Table 'project_member_role_use'
--
-- ---
DROP TABLE IF EXISTS project_member_role_use;
CREATE TABLE project_member_role_use (
	id serial NOT null PRIMARY KEY,
	project_role_id varchar(30) NULL,
	project_page_url_id int4 NOT NULL,
	project_role_search_yn bpchar(1) NULL,
	project_role_edit_yn bpchar(1) NULL,
	project_role_delete_yn bpchar(1) NULL,
	use_yn bpchar(1) NULL,
	creat_user varchar(30) NULL,
	creat_dt timestamp NULL,
	updt_user varchar(30) NULL,
	updt_dt timestamp NULL,
	project_role_modify_yn bpchar(1) NULL
);

-- ---
-- Table 'project_page_url'
--
-- ---
DROP TABLE IF EXISTS project_page_url;
CREATE TABLE project_page_url (
	id serial NOT null PRIMARY KEY,
	page_url varchar(255) NULL,
	page_url_name varchar(100) NULL,
	use_yn bpchar(1) NULL,
	creat_user varchar(30) NULL,
	creat_dt timestamp NULL,
	updt_user varchar(30) NULL,
	updt_dt timestamp NULL
);

-- ---
-- Table 'job_status'
--
-- ---
DROP TABLE IF EXISTS job_status;
CREATE TABLE job_status (
	id serial NOT null PRIMARY KEY,
	job_id int4 NOT NULL,
	status int4 NOT NULL,
	start_date timestamptz NULL,
	end_date timestamptz NULL
);

-- ---
-- Table 'dataset'
--
-- ---
DROP TABLE IF EXISTS dataset;
CREATE TABLE dataset (
	id serial NOT null PRIMARY KEY,
	member_id int4 NOT NULL,
	"name" varchar(100) NULL,
	description varchar(2000) NULL,
	use_yn bpchar(1) NULL,
	creat_user varchar(30) NULL,
	creat_dt timestamp NULL,
	updt_user varchar(30) NULL,
	updt_dt timestamp NULL
);

-- ---
-- Table 'dataset_detail'
--
-- ---
DROP TABLE IF EXISTS dataset_detail;
CREATE TABLE dataset_detail (
	id serial NOT null PRIMARY KEY,
	dataset_id int4 NOT NULL,
	molecule_id int4 NOT NULL,
	use_yn bpchar(1) NULL,
	creat_user varchar(30) NULL,
	creat_dt timestamp NULL,
	updt_user varchar(30) NULL,
	updt_dt timestamp NULL,
	add_dataset_id int4 NULL
);

-- ---
-- Table 'job_conformer'
--
-- ---
DROP TABLE IF EXISTS job_conformer;
CREATE TABLE job_conformer (
	id bigserial NOT null PRIMARY KEY,
  job_id BIGINT NOT NUll,
  conformer_id BIGINT NOT NULL,
  mol_3d VARCHAR(1000000) NOT NULL,
  "select" BOOLEAN DEFAULT false
);

-- ---
-- Table 'job_binding_site'
--
-- ---
DROP TABLE IF EXISTS job_binding_site;
CREATE TABLE job_binding_site (
	id bigserial NOT null PRIMARY KEY,
  job_id BIGINT NOT NUll,
  binding_site_id BIGINT NOT NULL
);

-- ---
-- Table 'job_conformer_cluster'
--
-- ---
DROP TABLE IF EXISTS job_conformer_cluster;
CREATE TABLE job_conformer_cluster (
	id bigserial NOT null PRIMARY KEY,
  job_conformer_id BIGINT NOT NUll,
  cluster_id BIGINT NOT NULL
);

-- ---
-- Table 'job_result'
--
-- ---
DROP TABLE IF EXISTS job_result;
CREATE TABLE job_result (
	id bigserial NOT null PRIMARY KEY,
  job_conformer_id BIGINT NOT NUll,
  index INTEGER,
  x varchar(20),
  value FLOAT,
  unit_id BIGINT NOT NULL
);

-- ---
-- Table 'unit'
--
-- ---
DROP TABLE IF EXISTS unit;
CREATE TABLE unit (
	id bigserial NOT null PRIMARY KEY,
  title VARCHAR NOT NUll
);

-- ---
-- Table 'pharmacophore'
--
-- ---
DROP TABLE IF EXISTS pharmacophore;
CREATE TABLE pharmacophore (
	id bigserial NOT null PRIMARY KEY,
  job_id BIGINT NOT NULL,
  file_path VARCHAR NOT NULL
);

-- ---
-- Table 'molecule_chembl'
--
-- ---
DROP TABLE IF EXISTS molecule_chembl;
CREATE TABLE molecule_chembl (
	id bigserial NOT null PRIMARY KEY,
  canonical_smiles molecule('sample')
);

-- ---
-- Table 'molecule_surechembl'
--
-- ---
DROP TABLE IF EXISTS molecule_surechembl;
CREATE TABLE molecule_surechembl (
	id bigserial NOT null PRIMARY KEY,
  canonical_smiles molecule('sample')
);

-- ---
-- Table 'fragment_chembl'
--
-- ---
DROP TABLE IF EXISTS fragment_chembl;
CREATE TABLE fragment_chembl (
	id bigserial NOT null PRIMARY KEY,
  canonical_smiles molecule('sample')
);

-- ---
-- Table 'fragment_surechembl'
--
-- ---
DROP TABLE IF EXISTS fragment_surechembl;
CREATE TABLE fragment_surechembl (
	id bigserial NOT null PRIMARY KEY,
  canonical_smiles molecule('sample')
);

-- ---
-- Table 'molecule_category'
--
-- ---
DROP TABLE IF EXISTS molecule_category;
CREATE TABLE molecule_category (
  id bigserial NOT NULL PRIMARY KEY,
  target_specificity INTEGER,
  patentability INTEGER
);

-- ---
-- Table 'dashboard'
--
-- ---
DROP TABLE IF EXISTS dashboard;
CREATE TABLE dashboard (
  id bigserial NOT NULL PRIMARY KEY,
  job_id BIGINT NOT NULL,
  log_date TIMESTAMP WITH TIME ZONE NOT NULL,
  file_path VARCHAR(1000)
);

-- ---
-- Index
-- ---

create index idx_compound_sourceid on compound(original_id);
create unique index idx_molecule_smiles on molecule(canonical_smiles);
create index idx_conformer_moleculeid on conformer(molecule_id);
create index idx_compmol_compoundid on compound_molecule(compound_id);
create index idx_compmol_moleculeid on compound_molecule(molecule_id);
create index idx_rgntcomp_compoundid on reagent_compound(compound_id);
create index idx_reagent_vendorid on reagent(vendor_id);
create index idx_reagent_code on reagent(code);
-- create index idx_pharmacophore_conformerid on pharmacophore(conformer_id);
--
-- create index idx_conformerfp_fp on conformer_fingerprint(fp);
-- create index idx_conformerfp_aa on conformer_fingerprint(aa);
-- create index idx_conformerfp_ad on conformer_fingerprint(ad);
-- create index idx_conformerfp_ah on conformer_fingerprint(ah);
-- create index idx_conformerfp_dd on conformer_fingerprint(dd);
-- create index idx_conformerfp_dh on conformer_fingerprint(dh);
-- create index idx_conformerfp_hh on conformer_fingerprint(hh);
--
-- create index idx_resultconf_jobid on result_conformer(job_id);
-- create index idx_resultmol_jobid on result_molecule(job_id);
--
create index idx_ligand_molid on ligand(molecule_id);
--create index idx_ligand_siteid on ligand(binding_site_id);
create index idx_site_ligandid on binding_site(ligand_id);
create index idx_site_proteinid on binding_site(protein_str_id);
create index idx_siteres_siteid on site_residue(binding_site_id);
create index idx_proteinstr_pdbid on protein_structure(pdb_id);
--
-- create index idx_datasetmol_datasetid on dataset_molecule(dataset_id);
--
create index idx_fragment_smiles on fragment(canonical_smiles);
create index idx_fragmentedge_start on fragment_edge(start);
create index idx_fragmentedge_end on fragment_edge("end");
create index idx_molfragment_molid on molecule_fragment(molecule_id);
create index idx_molfragment_fragid on molecule_fragment(fragment_id);
create index idx_scaffold_molfragid on scaffold(mol_frag_id);
--
create index idx_protein_code on protein_target(code);
create index idx_protein_pdb_proteinid on protein_pdb(protein_id);
create index idx_protein_pdb_strid on protein_pdb(protein_str_id);
create index idx_protein_name on protein_name(name);
create index idx_gene_name_name on gene_name(name);
create index idx_ppi_start on ppi(start);
create index idx_ppi_end on ppi("end");
create index idx_organism on organism(taxonomy_id);
create index idx_protein_organism_proteinid on protein_organism(protein_id);
create index idx_protein_organism_organismid on protein_organism(organism_id);
create index idx_disease_accession on disease(accession);
create index idx_disease_name on disease(name);
create index idx_protein_disease_proteinid on protein_disease(protein_id);
create index idx_protein_disease_diseaseid on protein_disease(disease_id);
create index idx_disease_synonym on disease_synonym(name);
create index idx_disease_reference_diseaseid on disease_reference(disease_id);
create index idx_disease_reference_referenceid on disease_reference(reference_id);
create index idx_ref_code on ref_code(code);
create index idx_keyword_name on keyword(name);
create index idx_keyword_go on keyword(go);
create index idx_protein_keyword_proteinid on protein_keyword(protein_id);
create index idx_protein_keyword_keywordid on protein_keyword(keyword_id);
create index idx_disease_keyword_diseaseid on disease_keyword(disease_id);
create index idx_disease_keyword_keywordid on disease_keyword(keyword_id);
create index idx_protein_uniprot_uniprotid on protein_uniprot(uniprot_id);
create index idx_protein_uniprot_proteinid on protein_uniprot(protein_id);
create index idx_protein_pdb_pdbid on pdb_chain(protein_pdb_id);
--
-- create index idx_resultassay_jobid on result_assay(job_id);
--
create index idx_physchem_mw on physchem(mw);
create index idx_physchem_logp on physchem(logp);
create index idx_physchem_tpsa on physchem(tpsa);
-- create index idx_lipinski_donors on lipinski(h_donors);
-- create index idx_lipinski_acceptors on lipinski(h_acceptors);
-- create index idx_lipinski_rbonds on lipinski(rot_bonds);
-- create index idx_lipinski_rings on lipinski(rings);
--
create index idx_conformer_id on atom3d(conformer_id);
create index idx_fp_all_01 on env_fingerprint(fp_all_01);
create index idx_fp_all_02 on env_fingerprint(fp_all_02);
create index idx_fp_all_03 on env_fingerprint(fp_all_03);
create index idx_fp_c_01 on env_fingerprint(fp_c_01);
create index idx_fp_c_02 on env_fingerprint(fp_c_02);
create index idx_fp_c_03 on env_fingerprint(fp_c_03);
create index idx_fp_n_01 on env_fingerprint(fp_n_01);
create index idx_fp_n_02 on env_fingerprint(fp_n_02);
create index idx_fp_n_03 on env_fingerprint(fp_n_03);
create index idx_fp_o_01 on env_fingerprint(fp_o_01);
create index idx_fp_o_02 on env_fingerprint(fp_o_02);
create index idx_fp_o_03 on env_fingerprint(fp_o_03);
create index idx_fp_p_01 on env_fingerprint(fp_p_01);
create index idx_fp_p_02 on env_fingerprint(fp_p_02);
create index idx_fp_p_03 on env_fingerprint(fp_p_03);
create index idx_fp_s_01 on env_fingerprint(fp_s_01);
create index idx_fp_s_02 on env_fingerprint(fp_s_02);
create index idx_fp_s_03 on env_fingerprint(fp_s_03);
create index idx_fp_hal_01 on env_fingerprint(fp_hal_01);
create index idx_fp_hal_02 on env_fingerprint(fp_hal_02);
create index idx_fp_hal_03 on env_fingerprint(fp_hal_03);
create index idx_feature_type on feature_fingerprint(type);
create index idx_feature_fp_triad on feature_fingerprint(fp_triad);
--
create index idx_gene_ensemblid on gene(ensembl_id);
create index idx_gene_name on gene(name);
create index idx_protein_target_geneid on protein_target(gene_id);
create index idx_hgnc_hgncid on hgnc(hgnc_id);
create index idx_hgnc_geneid on hgnc(gene_id);
--
create index idx_protein_target_chemblid on protein_target(chembl_id);
create index idx_target_molecule_proteinid on target_molecule(protein_id);
create index idx_target_molecule_moleculeid on target_molecule(molecule_id);
create index idx_target_molecule_activity on target_molecule(activity);
--
create index idx_protein_chembl_chemblid on protein_chembl(chembl_id);
--
create index idx_molecule_category_id on molecule_category(id);
--
create index idx_datasetdetail_datasetid on dataset_detail(dataset_id);
--
create index idx_dashboard_jobid on dashboard(job_id);

-- ---
-- Foreign Keys 
-- ---

ALTER TABLE compound ADD FOREIGN KEY (source_id) REFERENCES source (id);
ALTER TABLE compound_molecule ADD FOREIGN KEY (compound_id) REFERENCES compound (id);
ALTER TABLE compound_molecule ADD FOREIGN KEY (molecule_id) REFERENCES molecule (id);
ALTER TABLE reagent ADD FOREIGN KEY (vendor_id) REFERENCES vendor (id);
ALTER TABLE reagent_compound ADD FOREIGN KEY (compound_id) REFERENCES compound (id);
ALTER TABLE reagent_compound ADD FOREIGN KEY (id) REFERENCES reagent (id);
ALTER TABLE conformer ADD FOREIGN KEY (molecule_id) REFERENCES molecule (id);
-- ALTER TABLE pharmacophore ADD FOREIGN KEY (conformer_id) REFERENCES conformer (id);
-- ALTER TABLE pharmacophore ADD FOREIGN KEY (type_id) REFERENCES pharmacophore_type (id);
--
-- ALTER TABLE conformer_fingerprint ADD FOREIGN KEY (id) REFERENCES conformer (id);
--
--ALTER TABLE job_conformer ADD FOREIGN KEY (project_id) REFERENCES project (id);
--ALTER TABLE job_conformer ADD FOREIGN KEY (type) REFERENCES job_conformer_type (id);
-- ALTER TABLE result_conformer ADD FOREIGN KEY (job_id) REFERENCES job (id);
-- ALTER TABLE result_conformer ADD FOREIGN KEY (hit_id) REFERENCES conformer (id);
--
ALTER TABLE ligand ADD FOREIGN KEY (molecule_id) REFERENCES molecule (id);
ALTER TABLE binding_site ADD FOREIGN KEY (protein_str_id) REFERENCES protein_structure (id);
ALTER TABLE binding_site ADD FOREIGN KEY (ligand_id) REFERENCES ligand (id);
ALTER TABLE site_residue ADD FOREIGN KEY (residue_type_id) REFERENCES amino_acid (id);
ALTER TABLE site_residue ADD FOREIGN KEY (binding_site_id) REFERENCES binding_site (id);
--
ALTER TABLE molecule_fragment ADD FOREIGN KEY (molecule_id) REFERENCES molecule (id);
ALTER TABLE molecule_fragment ADD FOREIGN KEY (fragment_id) REFERENCES fragment (id);
ALTER TABLE fragment_edge ADD FOREIGN KEY (start) REFERENCES fragment (id);
ALTER TABLE fragment_edge ADD FOREIGN KEY ("end") REFERENCES fragment (id);
ALTER TABLE scaffold ADD FOREIGN KEY (mol_frag_id) REFERENCES molecule_fragment (id);
--ALTER TABLE job_molecule ADD FOREIGN KEY (project_id) REFERENCES project (id);
--ALTER TABLE job_molecule ADD FOREIGN KEY (type_id) REFERENCES job_molecule_type (id);
--ALTER TABLE result_molecule ADD FOREIGN KEY (job_id) REFERENCES job_molecule (id);
-- ALTER TABLE result_molecule ADD FOREIGN KEY (cluster_id) REFERENCES cluster (id);
-- ALTER TABLE result_molecule ADD FOREIGN KEY (molecule_id) REFERENCES molecule (id);
ALTER TABLE virtual ADD FOREIGN KEY (id) REFERENCES molecule (id);
--ALTER TABLE project ADD FOREIGN KEY (user_id) REFERENCES member (id);
--ALTER TABLE dataset ADD FOREIGN KEY (project_id) REFERENCES project (id);
-- ALTER TABLE dataset_molecule ADD FOREIGN KEY (dataset_id) REFERENCES dataset (id);
-- ALTER TABLE dataset_molecule ADD FOREIGN KEY (molecule_id) REFERENCES molecule (id);
ALTER TABLE cluster ADD FOREIGN KEY (fragment_id) REFERENCES fragment (id);
--
--ALTER TABLE custom ADD FOREIGN KEY (id) REFERENCES protein_structure (id);
--ALTER TABLE uniprot_code ADD FOREIGN KEY (protein_id) REFERENCES protein_target (id);
ALTER TABLE protein_name ADD FOREIGN KEY (protein_id) REFERENCES protein_target (id);
ALTER TABLE protein_pdb ADD FOREIGN KEY (protein_id) REFERENCES protein_target (id);
ALTER TABLE protein_pdb ADD FOREIGN KEY (protein_str_id) REFERENCES protein_structure (id);
ALTER TABLE protein_pdb ADD FOREIGN KEY (method) REFERENCES pdb_method (id);
ALTER TABLE gene_name ADD FOREIGN KEY (protein_id) REFERENCES protein_target (id);
ALTER TABLE protein_organism ADD FOREIGN KEY (protein_id) REFERENCES protein_target (id);
ALTER TABLE protein_organism ADD FOREIGN KEY (organism_id) REFERENCES organism (id);
ALTER TABLE protein_disease ADD FOREIGN KEY (protein_id) REFERENCES protein_target (id);
ALTER TABLE protein_disease ADD FOREIGN KEY (disease_id) REFERENCES disease (id);
ALTER TABLE keyword ADD FOREIGN KEY (category) REFERENCES keyword (id);
ALTER TABLE protein_keyword ADD FOREIGN KEY (protein_id) REFERENCES protein_target (id);
ALTER TABLE protein_keyword ADD FOREIGN KEY (keyword_id) REFERENCES keyword (id);
ALTER TABLE ppi ADD FOREIGN KEY (start) REFERENCES protein_target (id);
ALTER TABLE ppi ADD FOREIGN KEY ("end") REFERENCES protein_target (id);
ALTER TABLE ref_code ADD FOREIGN KEY (source_id) REFERENCES ref_source (id);
ALTER TABLE disease_reference ADD FOREIGN KEY (disease_id) REFERENCES disease (id);
ALTER TABLE disease_reference ADD FOREIGN KEY (reference_id) REFERENCES ref_code (id);
ALTER TABLE disease_synonym ADD FOREIGN KEY (disease_id) REFERENCES disease (id);
ALTER TABLE disease_keyword ADD FOREIGN KEY (disease_id) REFERENCES disease (id);
ALTER TABLE disease_keyword ADD FOREIGN KEY (keyword_id) REFERENCES keyword (id);
ALTER TABLE protein_uniprot ADD FOREIGN KEY (protein_id) REFERENCES protein_target (id);
ALTER TABLE protein_uniprot ADD FOREIGN KEY (uniprot_id) REFERENCES uniprot_code (id);
ALTER TABLE pdb_chain ADD FOREIGN KEY (protein_pdb_id) REFERENCES protein_pdb (id);
--
ALTER TABLE project_member ADD FOREIGN KEY (project_id) REFERENCES project (id);
ALTER TABLE project_member ADD FOREIGN KEY (member_id) REFERENCES member (id);
ALTER TABLE batch ADD FOREIGN KEY (project_id) REFERENCES project (id);
ALTER TABLE job ADD FOREIGN KEY (project_member_id) REFERENCES project_member (id);
ALTER TABLE job ADD FOREIGN KEY (batch_id) REFERENCES batch (id);
ALTER TABLE job ADD FOREIGN KEY (type_id) REFERENCES job_type (id);
ALTER TABLE job_edge ADD FOREIGN KEY (start) REFERENCES job (id);
ALTER TABLE job_edge ADD FOREIGN KEY ("end") REFERENCES job (id);
-- ALTER TABLE result_molecule ADD FOREIGN KEY (job_id) REFERENCES job (id);
-- ALTER TABLE result_conformer ADD FOREIGN KEY (job_id) REFERENCES job (id);
-- ALTER TABLE result_conformer ADD FOREIGN KEY (bsite_id) REFERENCES binding_site (id);
-- ALTER TABLE result_assay ADD FOREIGN KEY (job_id) REFERENCES job (id);
-- ALTER TABLE result_assay ADD FOREIGN KEY (molecule_id) REFERENCES molecule (id);
-- ALTER TABLE lead ADD FOREIGN KEY (result_assay_id) REFERENCES result_assay (id);
-- ALTER TABLE candidate ADD FOREIGN KEY (result_assay_id) REFERENCES result_assay (id);
--
ALTER TABLE nmap_model ADD FOREIGN KEY (disease_id) REFERENCES disease (id);
ALTER TABLE nmap_node ADD FOREIGN KEY (type_id) REFERENCES node_type (id);
ALTER TABLE nmap_edge ADD FOREIGN KEY (start) REFERENCES nmap_node (id);
ALTER TABLE nmap_edge ADD FOREIGN KEY ("end") REFERENCES nmap_node (id);
ALTER TABLE nmap_edge ADD FOREIGN KEY (model_id) REFERENCES nmap_model (id);
ALTER TABLE node_activity ADD FOREIGN KEY (network_id) REFERENCES nmap_network (id);
ALTER TABLE node_activity ADD FOREIGN KEY (node_id) REFERENCES nmap_node (id);
ALTER TABLE node_activity ADD FOREIGN KEY (activity_id) REFERENCES node_activity_type (id);
ALTER TABLE node_perturbation ADD FOREIGN KEY (node_id) REFERENCES nmap_node (id);
ALTER TABLE node_perturbation ADD FOREIGN KEY (activity_id) REFERENCES node_activity_type (id);
ALTER TABLE node_perturbation ADD FOREIGN KEY (perturbation_id) REFERENCES perturbation (id);
ALTER TABLE perturbation ADD FOREIGN KEY (network_id) REFERENCES nmap_network (id);
ALTER TABLE attractor ADD FOREIGN KEY (perturbation_id) REFERENCES perturbation (id);
ALTER TABLE attractor_phenotype ADD FOREIGN KEY (attractor_id) REFERENCES attractor (id);
ALTER TABLE attractor_phenotype ADD FOREIGN KEY (phenotype_id) REFERENCES nmap_phenotype (id);
ALTER TABLE network_cellline ADD FOREIGN KEY (network_id) REFERENCES nmap_network (id);
ALTER TABLE network_cellline ADD FOREIGN KEY (cellline_id) REFERENCES cellline (id);
--
-- ALTER TABLE target_node ADD FOREIGN KEY (protein_id) REFERENCES protein_target (id);
-- ALTER TABLE target_node ADD FOREIGN KEY (node_id) REFERENCES nmap_node (id);
--
ALTER TABLE tag_table ADD FOREIGN KEY (id) REFERENCES molecule (id);
--
--ALTER TABLE physchem ADD FOREIGN KEY (id) REFERENCES molecule (id);
--ALTER TABLE lipinski ADD FOREIGN KEY (id) REFERENCES molecule (id);
--ALTER TABLE admet ADD FOREIGN KEY (id) REFERENCES molecule (id);
--ALTER TABLE physchem ADD FOREIGN KEY (mw_class) REFERENCES mw_class_type (id);
--ALTER TABLE physchem ADD FOREIGN KEY (logp_class) REFERENCES logp_class_type (id);
--ALTER TABLE physchem ADD FOREIGN KEY (solubility_class) REFERENCES solubility_class_type (id);
--ALTER TABLE lipinski ADD FOREIGN KEY (h_donors_class) REFERENCES lipinski_class_type (id);
--ALTER TABLE lipinski ADD FOREIGN KEY (h_acceptors_class) REFERENCES lipinski_class_type (id);
--ALTER TABLE lipinski ADD FOREIGN KEY (rot_bonds_class) REFERENCES lipinski_class_type (id);
--ALTER TABLE lipinski ADD FOREIGN KEY (rings_class) REFERENCES lipinski_class_type (id);
--ALTER TABLE admet ADD FOREIGN KEY (caco2_permeability_class) REFERENCES permeability_class_type (id);
--ALTER TABLE admet ADD FOREIGN KEY (hia_class) REFERENCES hia_class_type (id);
--ALTER TABLE admet ADD FOREIGN KEY (metabolic_stability_class) REFERENCES stability_class_type (id);
--ALTER TABLE admet ADD FOREIGN KEY (ames_test_class) REFERENCES ames_class_type (id);
--ALTER TABLE admet ADD FOREIGN KEY (herg_inhibition_class) REFERENCES herg_class_type (id);
--
ALTER TABLE feature_fingerprint ADD FOREIGN KEY (conformer_id) REFERENCES conformer (id);
--
ALTER TABLE atom3d ADD FOREIGN KEY (conformer_id) REFERENCES conformer (id);
ALTER TABLE atom3d ADD FOREIGN KEY (element_id) REFERENCES element (id);
ALTER TABLE env_all ADD FOREIGN KEY (id) REFERENCES atom3d (id);
ALTER TABLE env_c ADD FOREIGN KEY (id) REFERENCES atom3d (id);
ALTER TABLE env_n ADD FOREIGN KEY (id) REFERENCES atom3d (id);
ALTER TABLE env_o ADD FOREIGN KEY (id) REFERENCES atom3d (id);
ALTER TABLE env_p ADD FOREIGN KEY (id) REFERENCES atom3d (id);
ALTER TABLE env_s ADD FOREIGN KEY (id) REFERENCES atom3d (id);
ALTER TABLE env_halogen ADD FOREIGN KEY (id) REFERENCES atom3d (id);
ALTER TABLE env_fingerprint ADD FOREIGN KEY (id) REFERENCES atom3d (id);
--
ALTER TABLE physchem ADD FOREIGN KEY (id) REFERENCES molecule (id);
ALTER TABLE absorption ADD FOREIGN KEY (id) REFERENCES molecule (id);
ALTER TABLE distribution ADD FOREIGN KEY (id) REFERENCES molecule (id);
ALTER TABLE metabolism ADD FOREIGN KEY (id) REFERENCES molecule (id);
ALTER TABLE toxicity ADD FOREIGN KEY (id) REFERENCES molecule (id);
ALTER TABLE leadlike_category ADD FOREIGN KEY (id) REFERENCES molecule (id);
ALTER TABLE leadlike_category ADD FOREIGN KEY (mw_class) REFERENCES mw_class_type (id);
ALTER TABLE leadlike_category ADD FOREIGN KEY (logp_class) REFERENCES logp_class_type (id);
ALTER TABLE leadlike_category ADD FOREIGN KEY (solubility_class) REFERENCES solubility_class_type (id);
ALTER TABLE leadlike_category ADD FOREIGN KEY (caco2_permeability_class) REFERENCES permeability_class_type (id);
ALTER TABLE leadlike_category ADD FOREIGN KEY (hia_class) REFERENCES hia_class_type (id);
ALTER TABLE leadlike_category ADD FOREIGN KEY (metabolic_stability_class) REFERENCES stability_class_type (id);
ALTER TABLE leadlike_category ADD FOREIGN KEY (ames_test_class) REFERENCES ames_class_type (id);
ALTER TABLE leadlike_category ADD FOREIGN KEY (herg_inhibition_class) REFERENCES herg_class_type (id);
--
ALTER TABLE hgnc ADD FOREIGN KEY (gene_id) REFERENCES gene (id);
ALTER TABLE gene_node ADD FOREIGN KEY (gene_id) REFERENCES gene (id);
ALTER TABLE gene_node ADD FOREIGN KEY (node_id) REFERENCES nmap_node (id);
ALTER TABLE protein_target ADD FOREIGN KEY (gene_id) REFERENCES gene (id);
--
ALTER TABLE target_molecule ADD FOREIGN KEY (protein_id) REFERENCES protein_target (id);
ALTER TABLE target_molecule ADD FOREIGN KEY (molecule_id) REFERENCES molecule (id);
--
ALTER TABLE protein_chembl ADD FOREIGN KEY (protein_id) REFERENCES protein_target(id);
-- 
ALTER TABLE job_binding_site ADD FOREIGN KEY (job_id) REFERENCES job (id);
ALTER TABLE job_binding_site ADD FOREIGN KEY (binding_site_id) REFERENCES binding_site (id);
ALTER TABLE job_conformer_cluster ADD FOREIGN KEY (job_conformer_id) REFERENCES job_conformer (id);
ALTER TABLE job_conformer_cluster ADD FOREIGN KEY (cluster_id) REFERENCES cluster (id);
ALTER TABLE job_conformer ADD FOREIGN KEY (conformer_id) REFERENCES conformer (id);
ALTER TABLE job_conformer ADD FOREIGN KEY (job_id) REFERENCES job (id);
ALTER TABLE job_result ADD FOREIGN KEY (job_conformer_id) REFERENCES job_conformer (id);
ALTER TABLE job_result ADD FOREIGN KEY (unit_id) REFERENCES unit (id);
ALTER TABLE pharmacophore ADD FOREIGN KEY (job_id) REFERENCES job (id);
ALTER TABLE lead ADD FOREIGN KEY (job_conformer_id) REFERENCES job_conformer (id);
ALTER TABLE candidate ADD FOREIGN KEY (job_conformer_id) REFERENCES job_conformer (id);
--
ALTER TABLE molecule_chembl ADD FOREIGN KEY (id) REFERENCES molecule (id);
ALTER TABLE molecule_surechembl ADD FOREIGN KEY (id) REFERENCES molecule (id);
ALTER TABLE fragment_chembl ADD FOREIGN KEY (id) REFERENCES fragment (id);
ALTER TABLE fragment_surechembl ADD FOREIGN KEY (id) REFERENCES fragment (id);
--
ALTER TABLE molecule_category ADD FOREIGN KEY (id) REFERENCES molecule (id);
--
ALTER TABLE dataset_detail ADD FOREIGN KEY (dataset_id) REFERENCES dataset (id);
ALTER TABLE dataset_detail ADD FOREIGN KEY (molecule_id) REFERENCES molecule (id);
ALTER TABLE job_status ADD FOREIGN KEY (job_id) REFERENCES job (id);
ALTER TABLE project_member_role_use ADD FOREIGN KEY (project_role_id) REFERENCES project_member_role (project_role_id);
ALTER TABLE project_member_role_use ADD FOREIGN KEY (project_page_url_id) REFERENCES project_page_url (id);
ALTER TABLE tc_menu ADD FOREIGN KEY (role_id) REFERENCES member_role (role_id);
ALTER TABLE tc_menu_mngr ADD FOREIGN KEY (menu_key) REFERENCES tc_menu (menu_key);
ALTER TABLE tc_menu_mngr ADD FOREIGN KEY (member_id) REFERENCES member (id);
ALTER TABLE tc_bbsctt ADD FOREIGN KEY (board_key) REFERENCES tc_bbs (board_key);
ALTER TABLE tc_bbsctt_cn ADD FOREIGN KEY (post_key) REFERENCES tc_bbsctt (post_key);
ALTER TABLE tc_bbsctt_atchmnfl ADD FOREIGN KEY (post_key) REFERENCES tc_bbsctt (post_key);
ALTER TABLE tc_bbsctt_comment ADD FOREIGN KEY (post_key) REFERENCES tc_bbsctt (post_key);
--
ALTER TABLE dashboard ADD FOREIGN KEY (job_id) REFERENCES job (id);

-- ---
-- Insert preset data
-- ---

-- source list
INSERT INTO source (name) VALUES ('ZINC');
INSERT INTO source (name) VALUES ('ChEMBL');
INSERT INTO source (name) VALUES ('SureChEMBL');
INSERT INTO source (name) VALUES ('ChemSpace');
INSERT INTO source (name) VALUES ('ChemDiv');
INSERT INTO source (name) VALUES ('Enamine');
INSERT INTO source (name) VALUES ('Chembridge');
INSERT INTO source (name) VALUES ('Lifechemicals');


-- pharmacophore_type
insert into pharmacophore_type (name) values ('HB ACCEPTOR');
insert into pharmacophore_type (name) values ('HB DONOR');
insert into pharmacophore_type (name) values ('HYDROPHOBIC');
insert into pharmacophore_type (name) values ('RING AROMATIC');
insert into pharmacophore_type (name) values ('POS IONIZABLE');
insert into pharmacophore_type (name) values ('NEG IONIZABLE');


-- vendor list
INSERT INTO vendor (name) VALUES ('ChemSpace');
INSERT INTO vendor (name) VALUES ('ChemDiv');
INSERT INTO vendor (name) VALUES ('Enamine');
INSERT INTO vendor (name) VALUES ('Chembridge');
INSERT INTO vendor (name) VALUES ('Lifechemicals');

-- distance_type
--insert into distance_type (name) values ('AA');
--insert into distance_type (name) values ('AD');
--insert into distance_type (name) values ('AH');
--insert into distance_type (name) values ('AN');
--insert into distance_type (name) values ('AP');
--insert into distance_type (name) values ('AR');
--insert into distance_type (name) values ('DD');
--insert into distance_type (name) values ('DH');
--insert into distance_type (name) values ('DN');
--insert into distance_type (name) values ('DP');
--insert into distance_type (name) values ('DR');
--insert into distance_type (name) values ('HH');
--insert into distance_type (name) values ('HN');
--insert into distance_type (name) values ('HP');
--insert into distance_type (name) values ('HR');
--insert into distance_type (name) values ('NN');
--insert into distance_type (name) values ('NP');
--insert into distance_type (name) values ('NR');
--insert into distance_type (name) values ('PP');
--insert into distance_type (name) values ('PR');
--insert into distance_type (name) values ('RR');

-- amino acid
insert into amino_acid (synonym, abbreviation, name) values ('', 'ABU', 'γ-Aminobutyric acid');
insert into amino_acid (synonym, abbreviation, name) values ('','ACD','Acidic unknown');
insert into amino_acid (synonym, abbreviation, name) values ('A', 'ALA', 'Alanine');
insert into amino_acid (synonym, abbreviation, name) values ('', 'ALB', 'β-Alanine');
insert into amino_acid (synonym, abbreviation, name) values ('', 'ALI', 'Aliphatic unknown');
insert into amino_acid (synonym, abbreviation, name) values ('R', 'ARG', 'Arginine');
insert into amino_acid (synonym, abbreviation, name) values ('', 'ARO', 'Aromatic unknown');
insert into amino_acid (synonym, abbreviation, name) values ('N', 'ASN', 'Asparagine');
insert into amino_acid (synonym, abbreviation, name) values ('D', 'ASP', 'Aspartic acid');
insert into amino_acid (synonym, abbreviation, name) values ('B', 'ASX', 'ASP/ASN ambiguous');
insert into amino_acid (synonym, abbreviation, name) values ('', 'BAS', 'Basic unknown');
insert into amino_acid (synonym, abbreviation, name) values ('C,CYH,CSH', 'CYS', 'Cysteine');
insert into amino_acid (synonym, abbreviation, name) values ('C,CSS,CYX', 'CYS', 'Cystine');
insert into amino_acid (synonym, abbreviation, name) values ('Q', 'GLN', 'Glutamine');
insert into amino_acid (synonym, abbreviation, name) values ('E', 'GLU', 'Glutamic acid');
insert into amino_acid (synonym, abbreviation, name) values ('Z', 'GLX', 'GLU/GLN ambiguous');
insert into amino_acid (synonym, abbreviation, name) values ('G', 'GLY', 'Glycine');
insert into amino_acid (synonym, abbreviation, name) values ('H', 'HIS', 'Histidine');
insert into amino_acid (synonym, abbreviation, name) values ('', 'HYP', 'Hydroxyproline');
insert into amino_acid (synonym, abbreviation, name) values ('I,ILU', 'ILE', 'Isoleucine');
insert into amino_acid (synonym, abbreviation, name) values ('L', 'LEU', 'Leucine');
insert into amino_acid (synonym, abbreviation, name) values ('K', 'LYS', 'Lysine');
insert into amino_acid (synonym, abbreviation, name) values ('M', 'MET', 'Methionine');
insert into amino_acid (synonym, abbreviation, name) values ('PGA', 'PCA', 'Pyrrolidone carboxylic acid');
insert into amino_acid (synonym, abbreviation, name) values ('F', 'PHE', 'Phenylalanine');
insert into amino_acid (synonym, abbreviation, name) values ('P,PRO,PRZ', 'PRO', 'Proline');
insert into amino_acid (synonym, abbreviation, name) values ('', 'SAR', 'Sarcosine');
insert into amino_acid (synonym, abbreviation, name) values ('S', 'SER', 'Serine');
insert into amino_acid (synonym, abbreviation, name) values ('T', 'THR', 'Threonine');
insert into amino_acid (synonym, abbreviation, name) values ('W,TRY', 'TRP', 'Tryptophan');
insert into amino_acid (synonym, abbreviation, name) values ('Y', 'TYR', 'Tryosine');
insert into amino_acid (synonym, abbreviation, name) values ('V', 'VAL', 'Valine');
insert into amino_acid (synonym, abbreviation, name) values ('', 'A', 'Adenosine');
insert into amino_acid (synonym, abbreviation, name) values ('', '1MA', '1-Methyladenosine');
insert into amino_acid (synonym, abbreviation, name) values ('', 'C', 'Cytidine');
insert into amino_acid (synonym, abbreviation, name) values ('', '5MC', '5-Methyladenosine');
insert into amino_acid (synonym, abbreviation, name) values ('', 'OMC', '2''-O-Methylcytidine');
insert into amino_acid (synonym, abbreviation, name) values ('', 'G', 'Guanosine');
insert into amino_acid (synonym, abbreviation, name) values ('', '1MG', '1-Methylguanosine');
insert into amino_acid (synonym, abbreviation, name) values ('', '2MG', 'N(2)-Methylguanosine');
insert into amino_acid (synonym, abbreviation, name) values ('', 'M2G', 'N(2)-Dimethylguanosine');
insert into amino_acid (synonym, abbreviation, name) values ('', '7MG', '7-Methylguanosine');
insert into amino_acid (synonym, abbreviation, name) values ('', 'OMG', '2''-O-Methylguanosine');
insert into amino_acid (synonym, abbreviation, name) values ('', 'YG', 'Wybutosine');
insert into amino_acid (synonym, abbreviation, name) values ('', 'I', 'Inosine');
insert into amino_acid (synonym, abbreviation, name) values ('', 'T', 'Thymidine');
insert into amino_acid (synonym, abbreviation, name) values ('', 'U', 'Uridine');
insert into amino_acid (synonym, abbreviation, name) values ('', '+U', 'Modified Uridine');
insert into amino_acid (synonym, abbreviation, name) values ('', 'H2U', 'Dihydrouridine');
insert into amino_acid (synonym, abbreviation, name) values ('', '5MU', 'Ribosylthymidine');
insert into amino_acid (synonym, abbreviation, name) values ('', 'PSU', 'Pseudouridine');
insert into amino_acid (synonym, abbreviation, name) values ('', 'ACE', 'Acetyl');
insert into amino_acid (synonym, abbreviation, name) values ('', 'FOR', 'Formyl');
insert into amino_acid (synonym, abbreviation, name) values ('H2O,WAT,OH2', 'HOH', 'Water');
insert into amino_acid (synonym, abbreviation, name) values ('', 'UNK', 'Unknown');
insert into amino_acid (synonym, abbreviation, name) values ('', 'DA', 'Deoxyadenosine');
insert into amino_acid (synonym, abbreviation, name) values ('', 'DC', 'Deoxycytidine');
insert into amino_acid (synonym, abbreviation, name) values ('', 'DG', 'Deoxyguanosine');
insert into amino_acid (synonym, abbreviation, name) values ('', 'DI', 'Deoxyinosine');
insert into amino_acid (synonym, abbreviation, name) values ('', 'DT', 'Deoxythymidine');
insert into amino_acid (synonym, abbreviation, name) values ('', 'DU', 'Deoxyuridine');
insert into amino_acid (synonym, abbreviation, name) values ('U', 'SEC', 'Selenocysteine');
insert into amino_acid (synonym, abbreviation, name) values ('O', 'PYL', 'Pyrrolysine');

--job type
insert into job_type(name) values ('screening1');
insert into job_type(name) values ('screening2');
insert into job_type(name) values ('docking');
insert into job_type(name) values ('clustering');
insert into job_type(name) values ('molecule_generation');
insert into job_type(name) values ('assay');

--node type
insert into node_type(name) values ('biological target');
insert into node_type(name) values ('stimulation');
insert into node_type(name) values ('phenotype');

--node activity type
insert into node_activity_type(name) values ('active');
insert into node_activity_type(name) values ('inactive');
insert into node_activity_type(name) values ('variable');

--nmap phenotype
insert into nmap_phenotype(name) values ('apoptosis');
insert into nmap_phenotype(name) values ('proliferation');

--nmap_model

--admet class
insert into logp_class_type(id, name) values (-2, 'Very lipophilic');
insert into logp_class_type(id, name) values (-1, 'Lipophilic');
insert into logp_class_type(id, name) values (0, 'Optimal');
insert into logp_class_type(id, name) values (1, 'Hydrophilic');
insert into logp_class_type(id, name) values (2, 'Very hydrophilic');
insert into logp_class_type(id, name) values (404, 'error');
insert into solubility_class_type(id, name) values (1, 'Soluble');
insert into solubility_class_type(id, name) values (0, 'Insoluble');
insert into solubility_class_type(id, name) values (-1, 'Highly insoluble');
insert into solubility_class_type(id, name) values (404, 'error');
insert into mw_class_type(id, name) values (1, 'Good');
insert into mw_class_type(id, name) values (0, 'Moderate');
insert into mw_class_type(id, name) values (-1, 'Bad');
insert into mw_class_type(id, name) values (404, 'error');
-- insert into lipinski_class_type(id, name) values (1, 'Good');
-- insert into lipinski_class_type(id, name) values (-1, 'Bad');
-- insert into lipinski_class_type(id, name) values (404, 'error');
insert into permeability_class_type(id, name) values (1, 'Highly permeable');
insert into permeability_class_type(id, name) values (0, 'Moderately permeable');
insert into permeability_class_type(id, name) values (-1, 'Poorly permeable');
insert into permeability_class_type(id, name) values (404, 'error');
insert into hia_class_type(id, name) values (1, 'Highly absorbed');
insert into hia_class_type(id, name) values (0, 'Moderately absorbed');
insert into hia_class_type(id, name) values (-1, 'Poorly absorbed');
insert into hia_class_type(id, name) values (404, 'error');
insert into stability_class_type(id, name) values (1, 'Stable in hlm');
insert into stability_class_type(id, name) values (0, 'Undefined');
insert into stability_class_type(id, name) values (-1, 'Unstable in HLM');
insert into stability_class_type(id, name) values (404, 'error');
insert into ames_class_type(id, name) values (1, 'Non-mutagenic');
insert into ames_class_type(id, name) values (0, 'Mutagenic');
insert into ames_class_type(id, name) values (-1, 'Nndefined');
insert into ames_class_type(id, name) values (404, 'error');
insert into herg_class_type(id, name) values (1, 'Non-inhibitor');
insert into herg_class_type(id, name) values (0, 'Inhibitor');
insert into herg_class_type(id, name) values (-1, 'Undefined');
insert into herg_class_type(id, name) values (404, 'error');

--element
INSERT INTO element (name, symbol, number) VALUES ('hydrogen', 'H', 1);
INSERT INTO element (name, symbol, number) VALUES ('boron', 'B', 5);
INSERT INTO element (name, symbol, number) VALUES ('carbon', 'C', 6);
INSERT INTO element (name, symbol, number) VALUES ('nitrogen', 'N', 7);
INSERT INTO element (name, symbol, number) VALUES ('oxygen', 'O', 8);
INSERT INTO element (name, symbol, number) VALUES ('fluorine', 'F', 9);
INSERT INTO element (name, symbol, number) VALUES ('sodium', 'Na', 11);
INSERT INTO element (name, symbol, number) VALUES ('magnesium', 'Mg', 12);
INSERT INTO element (name, symbol, number) VALUES ('silicon', 'Si', 14);
INSERT INTO element (name, symbol, number) VALUES ('phosphorus', 'P', 15);
INSERT INTO element (name, symbol, number) VALUES ('sulfur', 'S', 16);
INSERT INTO element  (name, symbol, number) VALUES ('chlorine', 'Cl', 17);
INSERT INTO element (name, symbol, number) VALUES ('potassium', 'K', 19);
INSERT INTO element (name, symbol, number) VALUES ('calcium', 'Ca', 20);
INSERT INTO element (name, symbol, number) VALUES ('arsenic', 'As',  33);
INSERT INTO element (name, symbol, number) VALUES ('bromine', 'Br', 35);
INSERT INTO element (name, symbol, number) VALUES ('tin', 'Sn', 50);
INSERT INTO element (name, symbol, number) VALUES ('iodine', 'I', 53);


-- ---
-- function
-- ---
CREATE OR REPLACE FUNCTION datediff(units character varying, start_t timestamp without time zone, end_t timestamp without time zone)
 RETURNS integer
 LANGUAGE plpgsql
AS $function$
   DECLARE
     diff_interval INTERVAL; 
     diff INT = 0;
     years_diff INT = 0;
   BEGIN
     IF units IN ('yy', 'yyyy', 'year', 'mm', 'm', 'month') THEN
       years_diff = DATE_PART('year', end_t) - DATE_PART('year', start_t);
 
       IF units IN ('yy', 'yyyy', 'year') THEN
         -- SQL Server does not count full years passed (only difference between year parts)
         RETURN years_diff;
       ELSE
         -- If end month is less than start month it will subtracted
         RETURN years_diff * 12 + (DATE_PART('month', end_t) - DATE_PART('month', start_t)); 
       END IF;
     END IF;
 
     -- Minus operator returns interval 'DDD days HH:MI:SS'  
     diff_interval = end_t - start_t;
 
     diff = diff + DATE_PART('day', diff_interval);
 
     IF units IN ('wk', 'ww', 'week') THEN
       diff = diff/7;
       RETURN diff;
     END IF;
 
     IF units IN ('dd', 'd', 'day') THEN
       RETURN diff;
     END IF;
 
     diff = diff * 24 + DATE_PART('hour', diff_interval); 
 
     IF units IN ('hh', 'hour') THEN
        RETURN diff;
     END IF;
 
     diff = diff * 60 + DATE_PART('minute', diff_interval);
 
     IF units IN ('mi', 'n', 'minute') THEN
        RETURN diff;
     END IF;
 
     diff = diff * 60 + DATE_PART('second', diff_interval);
 
     RETURN diff;
   END;
   $function$
;
