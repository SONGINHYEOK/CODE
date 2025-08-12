from django.contrib import admin
from .models import Molecule

@admin.register(Molecule)
class MoleculeAdmin(admin.ModelAdmin):
    list_display = ['id', 'smiles', 'molecular_weight', 'logp', 'source']
    search_fields = ['smiles', 'canonical_smiles', 'inchi_key']
    list_filter = ['source']
