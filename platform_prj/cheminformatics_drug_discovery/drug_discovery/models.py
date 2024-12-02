from django.db import models
from django.db.models import UniqueConstraint
from django.core.validators import MinValueValidator, MaxValueValidator
# from django.contrib.auth.models import User

class Target(models.Model):
    target_type = models.CharField(max_length=255, null=True, blank=True)
    description = models.CharField(max_length=255, null=True, blank=True)

    class Meta:
        db_table = "target"
    
    def __str__(self):
        return self.target_type

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
        return self.name

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
        return self.name

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
        return f"Job in Batch: {self.batch.name or 'N/A'}, Type: {self.job_type.name or 'N/A'}, Operator: {self.operator.member.name if self.operator else 'N/A'}"

class JobEdge(models.Model):
    start = models.ForeignKey(Job, on_delete=models.CASCADE, related_name="job_start_set", null=True, blank=True)
    end = models.ForeignKey(Job, on_delete=models.CASCADE, related_name="job_end_set", null=True, blank=True)

    class Meta:
        db_table = "job_edge"

    def __str__(self):
        return f"JobEdge from {self.start} to {self.end}"

class Bioassay(models.Model):
    name = models.CharField(max_length=255, null=True, blank=True)
    description = models.TextField(null=True, blank=True)
    protocol = models.TextField(null=True, blank=True)
    project = models.ForeignKey(Project, on_delete=models.CASCADE, null=True, blank=True)

    class Meta:
        db_table = "bioassay"
    
    def __str__(self):
        return self.name

class Gene(models.Model):
    gene_name = models.CharField(max_length=255, null=True, blank=True, unique=True)
    description = models.TextField(null=True, blank=True)

    class Meta:
        db_table = "gene"
    
    def __str__(self):
        return self.gene_name

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
            models.UniqueConstraint(fields=['gene', 'source_db', 'source_id'], name='unique_gene_source_constraint')
        ]

    def __str__(self):
        return f"GeneXref in Gene: {self.gene.gene_name}, source_id: {self.source_id}"


class Pathway(models.Model):
    internal_id = models.IntegerField(unique=True, blank=True)
    accession = models.CharField(max_length=255, unique=True, blank=True)

    class Meta:
        db_table = "pathway"
    
    def __str__(self):
        return f"internal_id: {self.internal_id}, accession: {self.accession}"

class TargetComment(models.Model):
    target = models.ForeignKey(Target, on_delete=models.CASCADE, null=True, blank=True)
    comment = models.TextField(null=True, blank=True)

    class Meta:
        db_table = "target_comment"
    
    def __str__(self):
        return self.comment

class TargetGene(models.Model):
    target = models.ForeignKey(Target, on_delete=models.CASCADE, null=True, blank=True)
    gene = models.ForeignKey(Gene, on_delete=models.CASCADE, null=True, blank=True)

    class Meta:
        db_table = "target_gene"
    
    def __str__(self):
        return self.gene.name

class TargetPathway(models.Model):
    target = models.ForeignKey(Target, on_delete=models.CASCADE, null=True, blank=True)
    pathway = models.ForeignKey(Pathway, on_delete=models.CASCADE, null=True, blank=True)

    class Meta:
        db_table = "target_pathway"
    
    def __str__(self):
        return f"TargetPathway in Target: {self.target.target_type}, Pathway: {self.pathway.id}"

class TargetPublication(models.Model):
    target = models.ForeignKey(Target, on_delete=models.CASCADE, null=True, blank=True)
    pubmed_id = models.IntegerField(null=True, blank=True)

    class Meta:
        db_table = "target_publication"
    
    def __str__(self):
        return self.pubmed_id

class Taxonomy(models.Model):
    ncbi_taxonomy = models.IntegerField(blank=True, unique=True)
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
        return self.ncbi_taxonomy

class ProteinXref(models.Model):
    ncbi_protein_accession = models.CharField(max_length=255, blank=True, unique=True)
    name = models.CharField(max_length=255, null=True, blank=True)
    ncbi_taxonomy = models.ForeignKey(Taxonomy, null=True, blank=True, on_delete=models.SET_NULL, to_field="ncbi_taxonomy")
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
        return self.name

class PDB(models.Model):
    original_id = models.CharField(max_length=10, unique=True)

    class Meta:
        db_table = "pdb"
    
    def __str__(self):
        return self.original_id

class Ligand(models.Model):
    pdb = models.ForeignKey(PDB, on_delete=models.SET_NULL, null=True, blank=True)

    class Meta:
        db_table = "ligand"
    
    def __str__(self):
        return f"Ligand in PDB: {self.pdb.original_id}"

class Molecule(models.Model):
    mol_3d = models.TextField(null=True, blank=True)

    class Meta:
        db_table = "molecule"

    def __str__(self):
        return self.mol_3d

class Compound(models.Model):
    ligand = models.ForeignKey(Ligand, on_delete=models.SET_NULL, null=True, blank=True)
    molecule = models.ForeignKey(Molecule, on_delete=models.SET_NULL, null=True, blank=True)
    molecular_formula = models.CharField(max_length=255, null=True, blank=True)
    molecular_weight = models.FloatField(null=True, blank=True)
    canonical_smiles = models.CharField(max_length=255, null=True, blank=True)
    isomeric_smiles = models.CharField(max_length=255, null=True, blank=True)

    class Meta:
        db_table = "compound"

    def __str__(self):
        return f"molecular_formula: {self.molecular_formula}"

class ProteinCompound(models.Model):
    protein = models.ForeignKey(Protein, on_delete=models.CASCADE, null=True, blank=True)
    compound = models.ForeignKey(Compound, on_delete=models.CASCADE, null=True, blank=True)

    class Meta:
        db_table = "protein_compound"

    def __str__(self):
        return f"ProteinCompound in Protein: {self.protein.name}, Compound: {self.compound.molecular_formula}"

class BioassayResult(models.Model):
    bioassay = models.ForeignKey(Bioassay, on_delete=models.CASCADE, null=True, blank=True)
    compound = models.ForeignKey(Compound, on_delete=models.CASCADE, null=True, blank=True)
    activity_outcome = models.CharField(max_length=255, null=True, blank=True)

    class Meta:
        db_table = "bioassay_result"
    
    def __str__(self):
        return self.activity_outcome

class CellLine(models.Model):
    internal_id = models.IntegerField(blank=True, unique=True)
    name = models.CharField(max_length=255, null=True, blank=True)
    tissue = models.CharField(max_length=255, null=True, blank=True)
    organism = models.CharField(max_length=255, null=True, blank=True)
    ncbi_taxonomy = models.ForeignKey(Taxonomy, null=True, blank=True, on_delete=models.SET_NULL, to_field="ncbi_taxonomy")
    synonyms = models.TextField(null=True, blank=True)

    class Meta:
        db_table = "cell_line"
        indexes = [
            models.Index(fields=['ncbi_taxonomy'], name='cell_line_ncbi_taxonomy_id_idx'),
        ]

    def __str__(self):
        return str(self.internal_id)

class Peptide(models.Model):
    sequence = models.CharField(max_length=255, null=True, blank=True)
    
    class Meta:
        db_table = "peptide"
    
    def __str__(self):
        return self.sequence

class ProteinPeptide(models.Model):
    protein = models.ForeignKey(Protein, on_delete=models.CASCADE, null=True, blank=True)
    peptide = models.ForeignKey(Peptide, on_delete=models.CASCADE, null=True, blank=True)

    class Meta:
        db_table = "protein_peptide"
    
    def __str__(self):
        return f"ProteinPeptide in Protein: {self.protein.name}, Peptide: {self.peptide.sequence}"

class JobCompound(models.Model):
    job = models.ForeignKey(Job, on_delete=models.CASCADE, null=True, blank=True)
    compound = models.ForeignKey(Compound, on_delete=models.CASCADE, null=True, blank=True)

    class Meta:
        db_table = "job_compound"
        constraints = [
            UniqueConstraint(fields=["job", "compound"], name="unique_job_compound")
        ]
    
    def __str__(self):
        return f"JobCompound in Job id: {self.job.id}, Job type: {self.job.job_type}, Compound: {self.compound.canonical_smiles}"

class JobPeptide(models.Model):
    job = models.ForeignKey(Job, on_delete=models.CASCADE, null=True, blank=True)
    peptide = models.ForeignKey(Peptide, on_delete=models.CASCADE, null=True, blank=True)

    class Meta:
        db_table = "job_peptide"
        constraints = [
            UniqueConstraint(fields=["job", "peptide"], name="unique_job_peptide")
        ]
    
    def __str__(self):
        return f"JobPeptide in Job id: {self.job.id}, Job type: {self.job.job_type}, Peptide: {self.peptide.sequence}" 

class CompoundValueType(models.Model):
    type_id = models.CharField(max_length=50)
    unit = models.CharField(max_length=50)

    class Meta:
        db_table = "compound_value_type"
    
    def __str__(self):
        return f"CompoundValueType type: {self.type_id}, unit: {self.unit}"

class PeptideValueType(models.Model):
    type_id = models.CharField(max_length=50)
    unit = models.CharField(max_length=50)

    class Meta:
        db_table = "peptide_value_type"
    
    def __str__(self):
        return f"PeptideValueType type: {self.type_id}, unit: {self.unit}"
    
class ResultCompound(models.Model):
    job_compound = models.ForeignKey(JobCompound, on_delete=models.CASCADE, null=True, blank=True)
    value = models.FloatField()
    value_type = models.ForeignKey(CompoundValueType, on_delete=models.CASCADE, null=True, blank=True)

    class Meta:
        db_table = "result_compound"
    
    def __str__(self):
        return f"ResultCompound in JobCompound job: {self.job_compound.job}, compound: {self.job_compound.compound}"

class ResultPeptide(models.Model):
    job = models.ForeignKey(Job, on_delete=models.CASCADE, null=True, blank=True)
    peptide = models.ForeignKey(Peptide, on_delete=models.CASCADE, null=True, blank=True)
    value = models.FloatField()
    value_type = models.ForeignKey(PeptideValueType, on_delete=models.CASCADE, null=True, blank=True)

    class Meta:
        db_table = "result_peptide"
        constraints = [
            UniqueConstraint(fields=["job", "peptide"], name="unique_result_peptide")
        ]
    
    def __str__(self):
        return f"ResultPeptide in Job id: {self.job.id}, Peptide: {self.peptide.sequence}"

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
        return f"Lipinski: {self.lipinski}, TPSA: {self.tpsa}"


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
        return f"HIA_Hou: {self.HIA_Hou}, Solubility_AqSolDB: {self.solubility_aqsol_db}"


class Distribution(models.Model):
    compound = models.OneToOneField(Compound, on_delete=models.CASCADE, primary_key=True)
    bbb_martins = models.SmallIntegerField(validators=[MinValueValidator(-1), MaxValueValidator(1)])
    ppbr_az = models.FloatField()
    vdss_lombardo = models.FloatField()

    class Meta:
        db_table = "distribution"

    def __str__(self):
        return f"BBB_Martins: {self.bbb_martins}, VDss_Lombardo: {self.vdss_lombardo}"


class Excretion(models.Model):
    compound = models.OneToOneField(Compound, on_delete=models.CASCADE, primary_key=True)
    half_life_obach = models.FloatField()
    clearance_hepatocyte_az = models.FloatField()
    clearance_microsome_az = models.FloatField()

    class Meta:
        db_table = "excretion"

    def __str__(self):
        return f"Half_Life_Obach: {self.half_life_obach}, Clearance_Hepatocyte_AZ: {self.clearance_hepatocyte_az}"


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
        return f"CYP1A2_Veith: {self.cyp1a2_veith}, CYP3A4_Veith: {self.cyp3a4_veith}"


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
        return f"hERG: {self.herg}, AMES: {self.ames}, DILI: {self.dili}"