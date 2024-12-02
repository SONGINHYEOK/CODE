from django.apps import AppConfig

class DrugDiscoveryConfig(AppConfig):
    default_auto_field = 'django.db.models.BigAutoField'
    name = 'drug_discovery'

    def ready(self):
        # signals import 제거
        pass