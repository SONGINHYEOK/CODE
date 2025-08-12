from django.contrib import admin
from .models import ScreeningJob, ScreeningResult

@admin.register(ScreeningJob)
class ScreeningJobAdmin(admin.ModelAdmin):
    list_display = ['id', 'batch', 'status', 'created_at', 'completed_at']
    list_filter = ['status', 'created_at']
    search_fields = ['batch__name']

@admin.register(ScreeningResult)
class ScreeningResultAdmin(admin.ModelAdmin):
    list_display = ['id', 'job', 'molecule', 'docking_score', 'is_hit']
    list_filter = ['is_hit', 'job']
    search_fields = ['molecule__smiles']
