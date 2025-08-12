from django.contrib import admin
from .models import Project, Target, Batch

@admin.register(Project)
class ProjectAdmin(admin.ModelAdmin):
    list_display = ['code', 'name', 'leader', 'status', 'created_at']
    list_filter = ['status', 'created_at']
    search_fields = ['code', 'name', 'description']

@admin.register(Target)
class TargetAdmin(admin.ModelAdmin):
    list_display = ['name', 'uniprot_id', 'pdb_id']
    search_fields = ['name', 'uniprot_id', 'pdb_id']

@admin.register(Batch)
class BatchAdmin(admin.ModelAdmin):
    list_display = ['name', 'project', 'target', 'status', 'created_at']
    list_filter = ['status', 'created_at']
    search_fields = ['name', 'project__name', 'target__name']
