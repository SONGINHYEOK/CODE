from django.db import models
from apps.core.models import TimeStampedModel
import uuid

class Molecule(TimeStampedModel):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    smiles = models.TextField()
    canonical_smiles = models.TextField(blank=True)
    inchi = models.TextField(blank=True)
    inchi_key = models.CharField(max_length=27, unique=True, null=True, blank=True)
    
    molecular_weight = models.FloatField(null=True, blank=True)
    logp = models.FloatField(null=True, blank=True)
    hbd = models.IntegerField(null=True, blank=True)
    hba = models.IntegerField(null=True, blank=True)
    
    source = models.CharField(max_length=50, blank=True)
    
    class Meta:
        db_table = 'molecules'
    
    def __str__(self):
        return self.smiles[:50]
