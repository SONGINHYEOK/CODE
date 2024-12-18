from django.db import models
from django.db.models import UniqueConstraint
from django.core.validators import MinValueValidator, MaxValueValidator

class Target(models.Model):
    target_type = models.CharField(max_length=255, null=True, blank=True)
    description = models.CharField(max_length=255, null=True, blank=True)

    class Meta:
        db_table = "target"
    
    def __str__(self):
        return self.target_type or "No Type"

class Project(models.Model):
    name = models.CharField(max_length=255)
    target = models.ForeignKey(Target, on_delete=models.CASCADE, null=True, blank=True)
    progress = models.CharField(max_length=255, null=True, blank=True)
    chem_or_pep = models.CharField(max_length=255, null=True, blank=True)

    class Meta:
        db_table = "project"
    
    def __str__(self):
        return self.name

class Member(models.Model):
    name = models.CharField(max_length=255, null=True, blank=True)

    class Meta:
        db_table = "member"

    def __str__(self):
        return self.name or "No Name"

class ProjectMember(models.Model):
    project = models.ForeignKey(Project, on_delete=models.CASCADE, null=True, blank=True)
    member = models.ForeignKey(Member, on_delete=models.CASCADE, null=True, blank=True)

    class Meta:
        db_table = "project_member"
        constraints = [
            UniqueConstraint(fields=['project', 'member'], name='unique_project_member')
        ]
    
    def __str__(self):
        return f"ProjectMember in Project: {self.project.id}, Member: {self.member.name}"

class Batch(models.Model):
    project = models.ForeignKey(Project, on_delete=models.CASCADE, null=True, blank=True)
    name = models.CharField(max_length=255, null=True, blank=True)

    class Meta:
        db_table = "batch"
    
    def __str__(self):
        return self.name or "No Name"

class JobType(models.Model):
    name = models.CharField(max_length=255, unique=True)

    class Meta:
        db_table = "job_type"
    
    def __str__(self):
        return self.name

class Job(models.Model):
    batch = models.ForeignKey(Batch, on_delete=models.CASCADE, null=True, blank=True)
    job_type = models.ForeignKey(JobType, on_delete=models.CASCADE, null=True, blank=True)
    operator = models.ForeignKey(ProjectMember, on_delete=models.CASCADE, null=True, blank=True)
    start_date = models.DateTimeField(null=True, blank=True)
    end_date = models.DateTimeField(null=True, blank=True)
    status = models.IntegerField(default=0)

    class Meta:
        db_table = "job"
    
    def __str__(self):
        return f"Job {self.id} - {self.job_type.name if self.job_type else 'No Type'}"

class JobEdge(models.Model):
    start = models.ForeignKey(Job, on_delete=models.CASCADE, related_name="job_start_set", null=True, blank=True)
    end = models.ForeignKey(Job, on_delete=models.CASCADE, related_name="job_end_set", null=True, blank=True)

    class Meta:
        db_table = "job_edge"

    def __str__(self):
        return f"JobEdge from {self.start.id} to {self.end.id}"

class Bioassay(models.Model):
    name = models.CharField(max_length=255, null=True, blank=True)
    description = models.TextField(null=True, blank=True)
    protocol = models.TextField(null=True, blank=True)
    project = models.ForeignKey(Project, on_delete=models.CASCADE, null=True, blank=True)

    class Meta:
        db_table = "bioassay"
    
    def __str__(self):
        return self.name or "No Name"

class Gene(models.Model):
    gene_name = models.CharField(max_length=255, unique=True, null=True, blank=True)
    description = models.TextField(null=True, blank=True)

    class Meta:
        db_table = "gene"
    
    def __str__(self):
        return self.gene_name or "No Name"

class GeneXref(models.Model):
    gene = models.ForeignKey(Gene, on_delete=models.CASCADE, null=True, blank=True)
    source_db = models.CharField(max_length=255, null=True, blank=True)
    source_id = models.CharField(max_length=255, null=True, blank=True)

    class Meta:
        db_table = "gene_xref"
        indexes = [
            models.Index(fields=['gene'])
        ]
        constraints = [
            models.UniqueConstraint(fields=['gene', 'source_db', 'source_id'], 
                                  name='unique_gene_source_constraint')
        ]

    def __str__(self):
        return f"GeneXref {self.source_db}:{self.source_id}"

class Pathway(models.Model):
    internal_id = models.IntegerField(unique=True, null=True, blank=True)
    accession = models.CharField(max_length=255, unique=True, null=True, blank=True)

    class Meta:
        db_table = "pathway"
    
    def __str__(self):
        return f"Pathway {self.accession}"

class TargetComment(models.Model):
    target = models.ForeignKey(Target, on_delete=models.CASCADE, null=True, blank=True)
    comment = models.TextField(null=True, blank=True)

    class Meta:
        db_table = "target_comment"
    
    def __str__(self):
        return self.comment or "No Comment"

class TargetGene(models.Model):
    target = models.ForeignKey(Target, on_delete=models.CASCADE, null=True, blank=True)
    gene = models.ForeignKey(Gene, on_delete=models.CASCADE, null=True, blank=True)

    class Meta:
        db_table = "target_gene"
    
    def __str__(self):
        return f"TargetGene {self.target.id}:{self.gene.gene_name}"

class TargetPathway(models.Model):
    target = models.ForeignKey(Target, on_delete=models.CASCADE, null=True, blank=True)
    pathway = models.ForeignKey(Pathway, on_delete=models.CASCADE, null=True, blank=True)

    class Meta:
        db_table = "target_pathway"
    
    def __str__(self):
        return f"TargetPathway {self.target.id}:{self.pathway.accession}"

class TargetPublication(models.Model):
    target = models.ForeignKey(Target, on_delete=models.CASCADE, null=True, blank=True)
    pubmed_id = models.IntegerField(null=True, blank=True)

    class Meta:
        db_table = "target_publication"
    
    def __str__(self):
        return str(self.pubmed_id)

class Taxonomy(models.Model):
    ncbi_taxonomy = models.IntegerField(unique=True, null=True, blank=True)
    scientific_name = models.CharField(max_length=255, null=True, blank=True)
    common_name = models.CharField(max_length=255, null=True, blank=True)
    synonyms = models.TextField(null=True, blank=True)
    domain = models.CharField(max_length=255, null=True, blank=True)
    rank = models.CharField(max_length=255, null=True, blank=True)
    rank_lineage = models.TextField(null=True, blank=True)
    parent_taxonomy = models.IntegerField(null=True, blank=True)

    class Meta:
        db_table = "taxonomy"
        indexes = [
            models.Index(fields=['ncbi_taxonomy'], name='idx_ncbi_taxonomy_id'),
        ]
    
    def __str__(self):
        return self.scientific_name or str(self.ncbi_taxonomy)

class PDB(models.Model):
    original_id = models.CharField(max_length=10, unique=True)

    class Meta:
        db_table = "pdb"
    
    def __str__(self):
        return self.original_id

class TargetPDB(models.Model):
    target = models.ForeignKey(Target, on_delete=models.CASCADE)
    pdb = models.ForeignKey(PDB, on_delete=models.CASCADE)

    class Meta:
        db_table = "target_pdb"
        constraints = [
            UniqueConstraint(fields=['target', 'pdb'], name='unique_target_pdb')
        ]

    def __str__(self):
        return f"TargetPDB {self.target.id}:{self.pdb.original_id}"

class Ligand(models.Model):
    pdb = models.ForeignKey(PDB, on_delete=models.SET_NULL, null=True, blank=True)

    class Meta:
        db_table = "ligand"
    
    def __str__(self):
        return f"Ligand {self.id}"

class JobPDB(models.Model):
    job = models.ForeignKey(Job, on_delete=models.CASCADE)
    pdb = models.ForeignKey(PDB, on_delete=models.CASCADE)

    class Meta:
        db_table = "job_pdb"
        constraints = [
            models.UniqueConstraint(
                fields=['job', 'pdb'],
                name='unique_job_pdb'
            )
        ]
    
    def __str__(self):
        return f"JobPDB in Job {self.job.id}, in PDB {self.pdb.id}"

class Molecule(models.Model):
    mol_3d = models.TextField(null=True, blank=True)

    class Meta:
        db_table = "molecule"

    def __str__(self):
        return f"Molecule {self.id}"

class Compound(models.Model):
    molecule = models.ForeignKey(Molecule, on_delete=models.SET_NULL, null=True, blank=True)
    molecular_formula = models.CharField(max_length=255, null=True, blank=True)
    molecular_weight = models.FloatField(null=True, blank=True)
    canonical_smiles = models.CharField(max_length=255, null=True, blank=True)
    isomeric_smiles = models.CharField(max_length=255, null=True, blank=True)

    class Meta:
        db_table = "compound"

    def __str__(self):
        return self.molecular_formula or f"Compound {self.id}"

class PDBLigandCompound(models.Model):
    ligand = models.ForeignKey(Ligand, on_delete=models.CASCADE)
    compound = models.ForeignKey(Compound, on_delete=models.CASCADE)

    class Meta:
        db_table = "pdbligand_compound"
        constraints = [
            UniqueConstraint(fields=['ligand', 'compound'], name='unique_ligand_compound')
        ]

    def __str__(self):
        return f"PDBLigandCompound {self.ligand.id}:{self.compound.id}"

class ProteinXref(models.Model):
    ncbi_protein_accession = models.CharField(max_length=255, unique=True, null=True, blank=True)
    name = models.CharField(max_length=255, null=True, blank=True)
    ncbi_taxonomy = models.ForeignKey(Taxonomy, null=True, blank=True, 
                                    on_delete=models.SET_NULL, to_field="ncbi_taxonomy")
    sequence = models.TextField(null=True, blank=True)
    synonyms = models.TextField(null=True, blank=True)

    class Meta:
        db_table = "protein_xref"
    
    def __str__(self):
        return self.ncbi_protein_accession or "No Accession"

class Protein(models.Model):
    name = models.CharField(max_length=255, null=True, blank=True)
    gene = models.ForeignKey(Gene, on_delete=models.SET_NULL, null=True, blank=True)
    organism = models.CharField(max_length=255, null=True, blank=True)
    sequence = models.TextField(null=True, blank=True)
    target = models.ForeignKey(Target, on_delete=models.CASCADE, null=True, blank=True)
    taxonomy = models.ForeignKey(Taxonomy, on_delete=models.CASCADE, null=True, blank=True)

    class Meta:
        db_table = "protein"
    
    def __str__(self):
        return self.name or f"Protein {self.id}"

class ProteinCompound(models.Model):
    protein = models.ForeignKey(Protein, on_delete=models.CASCADE, null=True, blank=True)
    compound = models.ForeignKey(Compound, on_delete=models.CASCADE, null=True, blank=True)

    class Meta:
        db_table = "protein_compound"

    def __str__(self):
        return f"ProteinCompound {self.protein.id}:{self.compound.id}"

class BioassayResult(models.Model):
    bioassay = models.ForeignKey(Bioassay, on_delete=models.CASCADE, null=True, blank=True)
    compound = models.ForeignKey(Compound, on_delete=models.CASCADE, null=True, blank=True)
    activity_outcome = models.CharField(max_length=255, null=True, blank=True)

    class Meta:
        db_table = "bioassay_result"
    
    def __str__(self):
        return f"BioassayResult {self.id}: {self.activity_outcome}"

class CellLine(models.Model):
    internal_id = models.IntegerField(unique=True, null=True, blank=True)
    name = models.CharField(max_length=255, null=True, blank=True)
    tissue = models.CharField(max_length=255, null=True, blank=True)
    organism = models.CharField(max_length=255, null=True, blank=True)
    ncbi_taxonomy = models.ForeignKey(Taxonomy, null=True, blank=True, 
                                    on_delete=models.SET_NULL, to_field="ncbi_taxonomy")
    synonyms = models.TextField(null=True, blank=True)

    class Meta:
        db_table = "cell_line"
        indexes = [
            models.Index(fields=['ncbi_taxonomy'], name='cell_line_ncbi_taxonomy_id_idx'),
        ]

    def __str__(self):
        return self.name or str(self.internal_id)

class Peptide(models.Model):
    sequence = models.CharField(max_length=255, null=True, blank=True)
    
    class Meta:
        db_table = "peptide"
    
    def __str__(self):
        return self.sequence or f"Peptide {self.id}"

class ProteinPeptide(models.Model):
    protein = models.ForeignKey(Protein, on_delete=models.CASCADE, null=True, blank=True)
    peptide = models.ForeignKey(Peptide, on_delete=models.CASCADE, null=True, blank=True)

    class Meta:
        db_table = "protein_peptide"
    
    def __str__(self):
        return f"ProteinPeptide {self.protein.id}:{self.peptide.id}"

class JobCompound(models.Model):
    job = models.ForeignKey(Job, on_delete=models.CASCADE, null=True, blank=True)
    compound = models.ForeignKey(Compound, on_delete=models.CASCADE, null=True, blank=True)

    class Meta:
        db_table = "job_compound"
        constraints = [
            UniqueConstraint(fields=["job", "compound"], name="unique_job_compound")
        ]
    
    def __str__(self):
        return f"JobCompound {self.job.id}:{self.compound.id}"

class JobPeptide(models.Model):
    job = models.ForeignKey(Job, on_delete=models.CASCADE, null=True, blank=True)
    peptide = models.ForeignKey(Peptide, on_delete=models.CASCADE, null=True, blank=True)

    class Meta:
        db_table = "job_peptide"
        constraints = [
            UniqueConstraint(fields=["job", "peptide"], name="unique_job_peptide")
        ]
    
    def __str__(self):
        return f"JobPeptide {self.job.id}:{self.peptide.id}"

class CompoundValueType(models.Model):
    type_id = models.CharField(max_length=50)
    unit = models.CharField(max_length=50)

    class Meta:
        db_table = "compound_value_type"
    
    def __str__(self):
        return f"CompoundValueType {self.type_id}:{self.unit}"

class PeptideValueType(models.Model):
    type_id = models.CharField(max_length=50)
    unit = models.CharField(max_length=50)

    class Meta:
        db_table = "peptide_value_type"
    
    def __str__(self):
        return f"PeptideValueType {self.type_id}:{self.unit}"
    
class ResultCompound(models.Model):
    job_compound = models.ForeignKey(JobCompound, on_delete=models.CASCADE, null=True, blank=True)
    value = models.FloatField()
    value_type = models.ForeignKey(CompoundValueType, on_delete=models.CASCADE, null=True, blank=True)

    class Meta:
        db_table = "result_compound"
    
    def __str__(self):
        return f"ResultCompound {self.job_compound.id}:{self.value}"

class ResultPeptide(models.Model):
    job = models.ForeignKey(Job, on_delete=models.CASCADE, null=True, blank=True)
    peptide = models.ForeignKey(Peptide, on_delete=models.CASCADE, null=True, blank=True)
    value = models.FloatField()
    value_type = models.ForeignKey(PeptideValueType, on_delete=models.CASCADE, null=True, blank=True)
    best_spot = models.IntegerField(null=True, blank=True)

    class Meta:
        db_table = "result_peptide"
        constraints = [
            UniqueConstraint(fields=["job", "peptide"], name="unique_result_peptide")
        ]
    
    def __str__(self):
        return f"ResultPeptide {self.job.id}:{self.peptide.id}:{self.value}"

class AssayType(models.Model):
    name = models.CharField(max_length=50)
    unit = models.CharField(max_length=20)
    description = models.TextField(null=True, blank=True)

    class Meta:
        db_table = "assay_type"
        constraints = [
            UniqueConstraint(fields=['name', 'unit'], name='unique_assay_type')
        ]

    def __str__(self):
        return f"{self.name} ({self.unit})"

class AssayResult(models.Model):
    job = models.ForeignKey(Job, on_delete=models.CASCADE)
    job_compound = models.ForeignKey(JobCompound, on_delete=models.CASCADE, null=True, blank=True)
    assay_type = models.ForeignKey(AssayType, on_delete=models.CASCADE, null=True, blank=True)
    value = models.FloatField()
    raw_data = models.JSONField(null=True, blank=True)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        db_table = "assay_result"
        constraints = [
            UniqueConstraint(
                fields=['job_compound', 'assay_type'], 
                name='unique_assay_result'
            )
        ]
        indexes = [
            models.Index(fields=['job'], name='idx_assay_result_job'),
            models.Index(
                fields=['job_compound', 'assay_type'], 
                name='idx_job_compound_assay_type'
            )
        ]

    def __str__(self):
        return f"AssayResult {self.job.id}:{self.assay_type.name}:{self.value}"

class AssayCondition(models.Model):
    assay_result = models.ForeignKey(AssayResult, on_delete=models.CASCADE)
    parameter = models.CharField(max_length=50)
    value = models.CharField(max_length=50)
    unit = models.CharField(max_length=20, null=True, blank=True)

    class Meta:
        db_table = "assay_condition"
        constraints = [
            UniqueConstraint(
                fields=['assay_result', 'parameter'], 
                name='unique_assay_condition'
            )
        ]

    def __str__(self):
        return f"{self.parameter}: {self.value}"

class Physicochemical(models.Model):
    compound = models.OneToOneField(Compound, on_delete=models.CASCADE, primary_key=True)
    molecular_weight = models.FloatField()
    logP = models.FloatField()
    hydrogen_bond_acceptors = models.IntegerField()
    hydrogen_bond_donors = models.IntegerField()
    lipinski = models.FloatField()
    qed = models.FloatField()
    stereo_centers = models.IntegerField()
    tpsa = models.FloatField()

    class Meta:
        db_table = "physicochemical"

    def __str__(self):
        return f"Physicochemical for Compound {self.compound.id}"

class Absorption(models.Model):
    compound = models.OneToOneField(Compound, on_delete=models.CASCADE, primary_key=True)
    HIA_Hou = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    bioavailability_ma = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    solubility_aqsol_db = models.FloatField()
    lipophilicity_astrazeneca = models.FloatField()
    hydration_free_energy_freesolv = models.FloatField()
    caco2_wang = models.FloatField()
    pampa_ncats = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    pgp_broccatelli = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])

    class Meta:
        db_table = "absorption"
    
    def __str__(self):
        return f"Absorption for Compound {self.compound.id}"

class Distribution(models.Model):
    compound = models.OneToOneField(Compound, on_delete=models.CASCADE, primary_key=True)
    bbb_martins = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    ppbr_az = models.FloatField()
    vdss_lombardo = models.FloatField()

    class Meta:
        db_table = "distribution"

    def __str__(self):
        return f"Distribution for Compound {self.compound.id}"

class Excretion(models.Model):
    compound = models.OneToOneField(Compound, on_delete=models.CASCADE, primary_key=True)
    half_life_obach = models.FloatField()
    clearance_hepatocyte_az = models.FloatField()
    clearance_microsome_az = models.FloatField()

    class Meta:
        db_table = "excretion"

    def __str__(self):
        return f"Excretion for Compound {self.compound.id}"

class Metabolism(models.Model):
    compound = models.OneToOneField(Compound, on_delete=models.CASCADE, primary_key=True)
    cyp1a2_veith = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    cyp2c19_veith = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    cyp2c9_veith = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    cyp2d6_veith = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    cyp3a4_veith = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    cyp2c9_substrate_carbonmangels = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    cyp2d6_substrate_carbonmangels = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    cyp3a4_substrate_carbonmangels = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])

    class Meta:
        db_table = "metabolism"

    def __str__(self):
        return f"Metabolism for Compound {self.compound.id}"
class ExperimentalAssayType(models.Model):
    name = models.CharField(max_length=50)  # 예: IC50, EC50, Ki 등
    unit = models.CharField(max_length=20)  # 예: nM, μM, % 등
    description = models.TextField(null=True, blank=True)

    class Meta:
        db_table = "experimental_assay_type"
        constraints = [
            models.UniqueConstraint(fields=['name', 'unit'], name='unique_experimental_assay_type')
        ]

    def __str__(self):
        return f"{self.name} ({self.unit})"
    
class ExperimentalAssayResult(models.Model):
    job = models.ForeignKey(Job, on_delete=models.CASCADE)
    job_compound = models.ForeignKey(JobCompound, on_delete=models.CASCADE, null=True, blank=True)
    assay_type = models.ForeignKey(ExperimentalAssayType, on_delete=models.CASCADE)
    value = models.FloatField()
    raw_data = models.JSONField(null=True, blank=True)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    
    class Meta:
        db_table = "experimental_assay_result"
        indexes = [
            models.Index(fields=['job'], name='idx_exp_assay_job'),
            models.Index(fields=['job_compound', 'assay_type'], 
                        name='idx_exp_job_compound_type')
        ]
        constraints = [
            models.UniqueConstraint(
                fields=['job_compound', 'assay_type'], 
                name='unique_experimental_assay_result'
            )
        ]
    
    def __str__(self):
        return f"ExperimentalAssayResult for Job {self.job.id} - {self.assay_type}: {self.value}"
    
class ExperimentalAssayCondition(models.Model):
    assay_result = models.ForeignKey(ExperimentalAssayResult, on_delete=models.CASCADE)
    parameter = models.CharField(max_length=50)  # 예: temperature, pH, time 등
    value = models.CharField(max_length=50)
    unit = models.CharField(max_length=20, null=True, blank=True)
    
    class Meta:
        db_table = "experimental_assay_condition"
        constraints = [
            models.UniqueConstraint(
                fields=['assay_result', 'parameter'], 
                name='unique_experimental_assay_condition'
            )
        ]

    def __str__(self):
        return f"{self.parameter}: {self.value} {self.unit or ''}"

class Toxicity(models.Model):
    compound = models.OneToOneField(Compound, on_delete=models.CASCADE, primary_key=True)
    herg = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    clintox = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    ames = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    dili = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    carcinogens_lagunin = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    ld50_zhu = models.FloatField()
    skin_reaction = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    nr_ar = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    nr_ar_lbd = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    nr_ahr = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    nr_aromatase = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    nr_er = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    nr_er_lbd = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    nr_ppar_gamma = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    sr_are = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    sr_atad5 = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    sr_hse = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    sr_mmp = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    sr_p53 = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])

    class Meta:
        db_table = "toxicity"

    def __str__(self):
        return f"Toxicity for Compound {self.compound.id}"

class Scaffold(models.Model):
    scaffold_smiles = models.CharField(max_length=255, null=True, blank=True)

    class Meta:
        db_table = "scaffold"
    
    def __str__(self):
        return f"Scaffold: {self.scaffold_smiles}"

class CompoundScaffold(models.Model):
    compound = models.ForeignKey(Compound, on_delete=models.CASCADE, null=True, blank=True)
    scaffold = models.ForeignKey(Scaffold, on_delete=models.CASCADE, null=True, blank=True)

    class Meta:
        db_table = "compound_scaffold"
        constraints = [
            UniqueConstraint(fields=['compound', 'scaffold'], name='unique_compound_scaffold')
        ]

    def __str__(self):
        return f"CompoundScaffold {self.compound.id}:{self.scaffold.id}"

class Dataset(models.Model):
    name = models.CharField(max_length=255, null=True, blank=True)
    owner = models.CharField(max_length=255, null=True, blank=True)
    count = models.IntegerField(null=True, blank=True)

    class Meta:
        db_table = "dataset"
    
    def __str__(self):
        return f"Dataset: {self.name}, {self.owner}, {self.count}"
        

class DatasetCompound(models.Model):
    dataset = models.ForeignKey(Dataset, on_delete=models.CASCADE)
    compound = models.ForeignKey(Compound, on_delete=models.CASCADE)

    class Meta:
        db_table = "dataset_compound"
        constraints = [
            UniqueConstraint(fields = ['dataset', 'compound'], name='unique_dataset_compound')
        ]
    
    def __str__(self):
        return f"DatasetCompound {self.dataset.id}: {self.compound.id}"

class DatasetCompoundDetail(models.Model):
    dataset_compound = models.ForeignKey(DatasetCompound, on_delete=models.CASCADE)
    source_result = models.ForeignKey(ResultCompound, on_delete=models.CASCADE)

    class Meta:
        db_table = "dataset_compound_detail"
        constraints = [
            models.UniqueConstraint(
                fields=['dataset_compound', 'source_result'],
                name='unique_dataset_compound_detail'
            )
        ]

    def __str__(self):
        return f"DatasetCompound - Detail: {self.dataset_compound_id}, Result: {self.source_result_id}"
