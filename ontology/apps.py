from django.apps import AppConfig
from django.db.models.signals import post_save


class OntologyConfig(AppConfig):
    name = 'ontology'

    def ready(self):
        from annotation.models import CachedWebResource
        from ontology.signals import gencc_post_save_handler

        post_save.connect(gencc_post_save_handler, sender=CachedWebResource)
