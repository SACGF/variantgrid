from django.apps import AppConfig
from django.db.models.signals import post_save


class OntologyConfig(AppConfig):
    name = 'ontology'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        # imported to activate receivers

        from annotation.models import CachedWebResource
        from ontology.models.ontology_search import (
            hpo_name_search,
            mondo_name_search,
            omim_name_search,
            ontology_search_hgnc,
            ontology_search_id,
        )

        # Registers receivers on import - noqa: F401 keeps the unused-import autofix from
        # silently unregistering them
        from ontology.signals import ontology_health_check, ontology_preview  # noqa: F401
        from ontology.signals.signals import gencc_post_save_handler
        from snpdb.search import search_registry
        # pylint: enable=import-outside-toplevel,unused-import

        search_registry.register(ontology_search_id, ontology_search_hgnc, omim_name_search, mondo_name_search,
                                 hpo_name_search)

        post_save.connect(gencc_post_save_handler, sender=CachedWebResource)
