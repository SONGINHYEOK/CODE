# Generated by Django 5.1 on 2024-09-29 03:17

import django.core.validators
import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Compound',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('molecular_formula', models.CharField(blank=True, max_length=255, null=True)),
                ('molecular_weight', models.FloatField(blank=True, null=True)),
                ('canonical_smiles', models.CharField(blank=True, max_length=255, null=True)),
                ('isomeric_smiles', models.CharField(blank=True, max_length=255, null=True)),
            ],
            options={
                'db_table': 'compound',
            },
        ),
        migrations.CreateModel(
            name='Batch',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(blank=True, max_length=255, null=True)),
            ],
            options={
                'db_table': 'batch',
            },
        ),
        migrations.CreateModel(
            name='Bioassay',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(blank=True, max_length=255, null=True)),
                ('description', models.TextField(blank=True, null=True)),
                ('protocol', models.TextField(blank=True, null=True)),
            ],
            options={
                'db_table': 'bioassay',
            },
        ),
        migrations.CreateModel(
            name='CompoundValueType',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('type_id', models.CharField(max_length=50)),
                ('unit', models.CharField(max_length=50)),
            ],
            options={
                'db_table': 'compound_value_type',
            },
        ),
        migrations.CreateModel(
            name='Gene',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('gene_name', models.CharField(blank=True, max_length=255, null=True, unique=True)),
                ('description', models.TextField(blank=True, null=True)),
            ],
            options={
                'db_table': 'gene',
            },
        ),
        migrations.CreateModel(
            name='JobType',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=255, unique=True)),
            ],
            options={
                'db_table': 'job_type',
            },
        ),
        migrations.CreateModel(
            name='Ligand',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
            ],
            options={
                'db_table': 'ligand',
            },
        ),
        migrations.CreateModel(
            name='Member',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(blank=True, max_length=255, null=True)),
            ],
            options={
                'db_table': 'member',
            },
        ),
        migrations.CreateModel(
            name='Molecule',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('mol_3d', models.TextField(blank=True, null=True)),
            ],
            options={
                'db_table': 'molecule',
            },
        ),
        migrations.CreateModel(
            name='Pathway',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('internal_id', models.IntegerField(blank=True, unique=True)),
                ('accession', models.CharField(blank=True, max_length=255, unique=True)),
            ],
            options={
                'db_table': 'pathway',
            },
        ),
        migrations.CreateModel(
            name='PDB',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('original_id', models.CharField(max_length=10, unique=True)),
            ],
            options={
                'db_table': 'pdb',
            },
        ),
        migrations.CreateModel(
            name='Peptide',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('sequence', models.CharField(blank=True, max_length=255, null=True)),
            ],
            options={
                'db_table': 'peptide',
            },
        ),
        migrations.CreateModel(
            name='PeptideValueType',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('type_id', models.CharField(max_length=50)),
                ('unit', models.CharField(max_length=50)),
            ],
            options={
                'db_table': 'peptide_value_type',
            },
        ),
        migrations.CreateModel(
            name='Project',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=255)),
                ('progress', models.CharField(blank=True, max_length=255, null=True)),
                ('chem_or_pep', models.CharField(blank=True, max_length=255, null=True)),
            ],
            options={
                'db_table': 'project',
            },
        ),
        migrations.CreateModel(
            name='Target',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('target_type', models.CharField(blank=True, max_length=255, null=True)),
                ('description', models.CharField(blank=True, max_length=255, null=True)),
            ],
            options={
                'db_table': 'target',
            },
        ),
        migrations.CreateModel(
            name='Absorption',
            fields=[
                ('compound', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, primary_key=True, serialize=False, to='drug_discovery.compound')),
                ('HIA_Hou', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('bioavailability_ma', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('solubility_aqsol_db', models.FloatField()),
                ('lipophilicity_astrazeneca', models.FloatField()),
                ('hydration_free_energy_freesolv', models.FloatField()),
                ('caco2_wang', models.FloatField()),
                ('pampa_ncats', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('pgp_broccatelli', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
            ],
            options={
                'db_table': 'absorption',
            },
        ),
        migrations.CreateModel(
            name='Distribution',
            fields=[
                ('compound', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, primary_key=True, serialize=False, to='drug_discovery.compound')),
                ('bbb_martins', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('ppbr_az', models.FloatField()),
                ('vdss_lombardo', models.FloatField()),
            ],
            options={
                'db_table': 'distribution',
            },
        ),
        migrations.CreateModel(
            name='Excretion',
            fields=[
                ('compound', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, primary_key=True, serialize=False, to='drug_discovery.compound')),
                ('half_life_obach', models.FloatField()),
                ('clearance_hepatocyte_az', models.FloatField()),
                ('clearance_microsome_az', models.FloatField()),
            ],
            options={
                'db_table': 'excretion',
            },
        ),
        migrations.CreateModel(
            name='Metabolism',
            fields=[
                ('compound', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, primary_key=True, serialize=False, to='drug_discovery.compound')),
                ('cyp1a2_veith', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('cyp2c19_veith', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('cyp2c9_veith', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('cyp2d6_veith', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('cyp3a4_veith', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('cyp2c9_substrate_carbonmangels', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('cyp2d6_substrate_carbonmangels', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('cyp3a4_substrate_carbonmangels', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
            ],
            options={
                'db_table': 'metabolism',
            },
        ),
        migrations.CreateModel(
            name='Physicochemical',
            fields=[
                ('compound', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, primary_key=True, serialize=False, to='drug_discovery.compound')),
                ('molecular_weight', models.FloatField()),
                ('logP', models.FloatField()),
                ('hydrogen_bond_acceptors', models.IntegerField()),
                ('hydrogen_bond_donors', models.IntegerField()),
                ('lipinski', models.FloatField()),
                ('qed', models.FloatField()),
                ('stereo_centers', models.IntegerField()),
                ('tpsa', models.FloatField()),
            ],
            options={
                'db_table': 'physicochemical',
            },
        ),
        migrations.CreateModel(
            name='Toxicity',
            fields=[
                ('compound', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, primary_key=True, serialize=False, to='drug_discovery.compound')),
                ('herg', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('clintox', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('ames', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('dili', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('carcinogens_lagunin', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('ld50_zhu', models.FloatField()),
                ('skin_reaction', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('nr_ar', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('nr_ar_lbd', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('nr_ahr', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('nr_aromatase', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('nr_er', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('nr_er_lbd', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('nr_ppar_gamma', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('sr_are', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('sr_atad5', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('sr_hse', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('sr_mmp', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
                ('sr_p53', models.SmallIntegerField(validators=[django.core.validators.MinValueValidator(-1), django.core.validators.MaxValueValidator(1)])),
            ],
            options={
                'db_table': 'toxicity',
            },
        ),
        migrations.CreateModel(
            name='BioassayResult',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('activity_outcome', models.CharField(blank=True, max_length=255, null=True)),
                ('bioassay', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.bioassay')),
                ('compound', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.compound')),
            ],
            options={
                'db_table': 'bioassay_result',
            },
        ),
        migrations.CreateModel(
            name='Job',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('start_date', models.DateTimeField(blank=True, null=True)),
                ('end_date', models.DateTimeField(blank=True, null=True)),
                ('status', models.IntegerField(default=0)),
                ('batch', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.batch')),
                ('job_type', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.jobtype')),
            ],
            options={
                'db_table': 'job',
            },
        ),
        migrations.CreateModel(
            name='JobCompound',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('compound', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.compound')),
                ('job', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.job')),
            ],
            options={
                'db_table': 'job_compound',
            },
        ),
        migrations.CreateModel(
            name='JobEdge',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('end', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='job_end_set', to='drug_discovery.job')),
                ('start', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='job_start_set', to='drug_discovery.job')),
            ],
            options={
                'db_table': 'job_edge',
            },
        ),
        migrations.AddField(
            model_name='compound',
            name='ligand',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, to='drug_discovery.ligand'),
        ),
        migrations.AddField(
            model_name='compound',
            name='molecule',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, to='drug_discovery.molecule'),
        ),
        migrations.AddField(
            model_name='ligand',
            name='pdb',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, to='drug_discovery.pdb'),
        ),
        migrations.CreateModel(
            name='JobPeptide',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('job', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.job')),
                ('peptide', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.peptide')),
            ],
            options={
                'db_table': 'job_peptide',
            },
        ),
        migrations.AddField(
            model_name='bioassay',
            name='project',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.project'),
        ),
        migrations.AddField(
            model_name='batch',
            name='project',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.project'),
        ),
        migrations.CreateModel(
            name='ProjectMember',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('member', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.member')),
                ('project', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.project')),
            ],
            options={
                'db_table': 'project_member',
            },
        ),
        migrations.AddField(
            model_name='job',
            name='operator',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.projectmember'),
        ),
        migrations.CreateModel(
            name='Protein',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(blank=True, max_length=255, null=True)),
                ('organism', models.CharField(blank=True, max_length=255, null=True)),
                ('sequence', models.TextField(blank=True, null=True)),
                ('gene', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, to='drug_discovery.gene')),
                ('target', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.target')),
            ],
            options={
                'db_table': 'protein',
            },
        ),
        migrations.CreateModel(
            name='ProteinCompound',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('compound', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.compound')),
                ('protein', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.protein')),
            ],
            options={
                'db_table': 'protein_compound',
            },
        ),
        migrations.CreateModel(
            name='ProteinPeptide',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('peptide', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.peptide')),
                ('protein', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.protein')),
            ],
            options={
                'db_table': 'protein_peptide',
            },
        ),
        migrations.CreateModel(
            name='ResultCompound',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('value', models.FloatField()),
                ('job_compound', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.jobcompound')),
                ('value_type', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.compoundvaluetype')),
            ],
            options={
                'db_table': 'result_compound',
            },
        ),
        migrations.CreateModel(
            name='ResultPeptide',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('value', models.FloatField()),
                ('job', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.job')),
                ('peptide', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.peptide')),
                ('value_type', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.peptidevaluetype')),
            ],
            options={
                'db_table': 'result_peptide',
            },
        ),
        migrations.AddField(
            model_name='project',
            name='target',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.target'),
        ),
        migrations.CreateModel(
            name='TargetComment',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('comment', models.TextField(blank=True, null=True)),
                ('target', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.target')),
            ],
            options={
                'db_table': 'target_comment',
            },
        ),
        migrations.CreateModel(
            name='TargetGene',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('gene', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.gene')),
                ('target', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.target')),
            ],
            options={
                'db_table': 'target_gene',
            },
        ),
        migrations.CreateModel(
            name='TargetPathway',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('pathway', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.pathway')),
                ('target', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.target')),
            ],
            options={
                'db_table': 'target_pathway',
            },
        ),
        migrations.CreateModel(
            name='TargetPublication',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('pubmed_id', models.IntegerField(blank=True, null=True)),
                ('target', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.target')),
            ],
            options={
                'db_table': 'target_publication',
            },
        ),
        migrations.CreateModel(
            name='Taxonomy',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('ncbi_taxonomy', models.IntegerField(blank=True, unique=True)),
                ('scientific_name', models.CharField(blank=True, max_length=255, null=True)),
                ('common_name', models.CharField(blank=True, max_length=255, null=True)),
                ('synonyms', models.TextField(blank=True, null=True)),
                ('domain', models.CharField(blank=True, max_length=255, null=True)),
                ('rank', models.CharField(blank=True, max_length=255, null=True)),
                ('rank_lineage', models.TextField(blank=True, null=True)),
                ('parent_taxonomy', models.IntegerField(blank=True, null=True)),
            ],
            options={
                'db_table': 'taxonomy',
                'indexes': [models.Index(fields=['ncbi_taxonomy'], name='idx_ncbi_taxonomy_id')],
            },
        ),
        migrations.CreateModel(
            name='ProteinXref',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('ncbi_protein_accession', models.CharField(blank=True, max_length=255, unique=True)),
                ('name', models.CharField(blank=True, max_length=255, null=True)),
                ('sequence', models.TextField(blank=True, null=True)),
                ('synonyms', models.TextField(blank=True, null=True)),
                ('ncbi_taxonomy', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, to='drug_discovery.taxonomy', to_field='ncbi_taxonomy')),
            ],
            options={
                'db_table': 'protein_xref',
            },
        ),
        migrations.AddField(
            model_name='protein',
            name='taxonomy',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.taxonomy'),
        ),
        migrations.CreateModel(
            name='CellLine',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('internal_id', models.IntegerField(blank=True, unique=True)),
                ('name', models.CharField(blank=True, max_length=255, null=True)),
                ('tissue', models.CharField(blank=True, max_length=255, null=True)),
                ('organism', models.CharField(blank=True, max_length=255, null=True)),
                ('synonyms', models.TextField(blank=True, null=True)),
                ('ncbi_taxonomy', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, to='drug_discovery.taxonomy', to_field='ncbi_taxonomy')),
            ],
            options={
                'db_table': 'cell_line',
            },
        ),
        migrations.CreateModel(
            name='GeneXref',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('source_db', models.CharField(blank=True, max_length=255, null=True)),
                ('source_id', models.CharField(blank=True, max_length=255, null=True)),
                ('gene', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='drug_discovery.gene')),
            ],
            options={
                'db_table': 'gene_xref',
                'indexes': [models.Index(fields=['gene'], name='gene_xref_gene_id_6bfb14_idx')],
                'constraints': [models.UniqueConstraint(fields=('gene', 'source_db', 'source_id'), name='unique_gene_source_constraint')],
            },
        ),
        migrations.AddConstraint(
            model_name='jobcompound',
            constraint=models.UniqueConstraint(fields=('job', 'compound'), name='unique_job_compound'),
        ),
        migrations.AddConstraint(
            model_name='jobpeptide',
            constraint=models.UniqueConstraint(fields=('job', 'peptide'), name='unique_job_peptide'),
        ),
        migrations.AddConstraint(
            model_name='projectmember',
            constraint=models.UniqueConstraint(fields=('project', 'member'), name='unique_project_member'),
        ),
        migrations.AddConstraint(
            model_name='resultpeptide',
            constraint=models.UniqueConstraint(fields=('job', 'peptide'), name='unique_result_peptide'),
        ),
        migrations.AddIndex(
            model_name='cellline',
            index=models.Index(fields=['ncbi_taxonomy'], name='cell_line_ncbi_taxonomy_id_idx'),
        ),
    ]