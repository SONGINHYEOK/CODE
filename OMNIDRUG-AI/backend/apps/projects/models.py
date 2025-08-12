from django.db import models
from django.contrib.auth import get_user_model
from apps.core.models import TimeStampedModel
import uuid

User = get_user_model()

class Project(TimeStampedModel):
    STATUS_CHOICES = [
        ('planning', 'Planning'),
        ('active', 'Active'),
        ('completed', 'Completed'),
        ('archived', 'Archived'),
    ]
    
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    code = models.CharField(max_length=50, unique=True)
    name = models.CharField(max_length=255)
    description = models.TextField(blank=True)
    leader = models.ForeignKey(User, on_delete=models.PROTECT, related_name='led_projects')
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='planning')
    
    therapeutic_area = models.CharField(max_length=100, blank=True)
    target_indication = models.CharField(max_length=255, blank=True)
    
    class Meta:
        db_table = 'projects'
        ordering = ['-created_at']
    
    def __str__(self):
        return f"{self.code}: {self.name}"

class Target(TimeStampedModel):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=255)
    uniprot_id = models.CharField(max_length=20, unique=True, null=True, blank=True)
    pdb_id = models.CharField(max_length=10, blank=True)
    sequence = models.TextField(blank=True)
    
    class Meta:
        db_table = 'targets'
    
    def __str__(self):
        return self.name

class Batch(TimeStampedModel):
    STATUS_CHOICES = [
        ('pending', 'Pending'),
        ('running', 'Running'),
        ('completed', 'Completed'),
        ('failed', 'Failed'),
    ]
    
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    project = models.ForeignKey(Project, on_delete=models.CASCADE, related_name='batches')
    name = models.CharField(max_length=255)
    target = models.ForeignKey(Target, on_delete=models.PROTECT)
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='pending')
    
    # UniDock parameters
    center_x = models.FloatField(null=True, blank=True)
    center_y = models.FloatField(null=True, blank=True)
    center_z = models.FloatField(null=True, blank=True)
    size_x = models.FloatField(default=20.0)
    size_y = models.FloatField(default=20.0)
    size_z = models.FloatField(default=20.0)
    
    class Meta:
        db_table = 'batches'
        unique_together = ['project', 'name']
    
    def __str__(self):
        return f"{self.project.code} - {self.name}"
