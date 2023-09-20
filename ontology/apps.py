from django.apps import AppConfig
from django.db.models.signals import post_save


class OntologyConfig(AppConfig):
    name = 'ontology'

    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        # imported to activate receivers
        from ontology.signals import ontology_health_check
        from ontology.signals import ontology_preview

        from annotation.models import CachedWebResource
        from ontology.signals.signals import gencc_post_save_handler
        # pylint: enable=import-outside-toplevel,unused-import

        post_save.connect(gencc_post_save_handler, sender=CachedWebResource)
