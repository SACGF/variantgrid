from django.apps import AppConfig
from django.db.models.signals import post_save


class OntologyConfig(AppConfig):
    name = 'ontology'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        # imported to activate receivers

        from annotation.models import CachedWebResource
        from ontology.signals.signals import gencc_post_save_handler
        # pylint: enable=import-outside-toplevel,unused-import

        post_save.connect(gencc_post_save_handler, sender=CachedWebResource)
