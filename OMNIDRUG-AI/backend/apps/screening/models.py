from django.db import models
from apps.core.models import TimeStampedModel
from apps.projects.models import Batch
from apps.molecules.models import Molecule
import uuid

class ScreeningJob(TimeStampedModel):
    STATUS_CHOICES = [
        ('queued', 'Queued'),
        ('running', 'Running'),
        ('completed', 'Completed'),
        ('failed', 'Failed'),
    ]
    
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    batch = models.ForeignKey(Batch, on_delete=models.CASCADE, related_name='screening_jobs')
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='queued')
    
    # Input files
    receptor_file = models.FileField(upload_to='receptors/', null=True, blank=True)
    ligand_file = models.FileField(upload_to='ligands/', null=True, blank=True)
    
    # UniDock parameters
    exhaustiveness = models.IntegerField(default=8)
    num_modes = models.IntegerField(default=9)
    energy_range = models.FloatField(default=3.0)
    
    # Results
    result_file = models.FileField(upload_to='results/', null=True, blank=True)
    
    # Task tracking
    celery_task_id = models.CharField(max_length=255, blank=True)
    started_at = models.DateTimeField(null=True, blank=True)
    completed_at = models.DateTimeField(null=True, blank=True)
    error_message = models.TextField(blank=True)
    
    class Meta:
        db_table = 'screening_jobs'

class ScreeningResult(TimeStampedModel):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    job = models.ForeignKey(ScreeningJob, on_delete=models.CASCADE, related_name='results')
    molecule = models.ForeignKey(Molecule, on_delete=models.CASCADE)
    
    docking_score = models.FloatField()
    rmsd = models.FloatField(null=True, blank=True)
    binding_pose = models.JSONField(default=dict)
    
    is_hit = models.BooleanField(default=False)
    
    class Meta:
        db_table = 'screening_results'
        ordering = ['docking_score']
