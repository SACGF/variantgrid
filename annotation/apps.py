from django.apps.config import AppConfig
from django.db.models.signals import post_save


class AnnotationConfig(AppConfig):
    name = 'annotation'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        # imported to activate receivers

        from Bio import Entrez
        from django.conf import settings
        from annotation.models import CachedWebResource
        from annotation.signals.manual_signals import clinvar_citations_post_save_handler
        from annotation.signals import citation_preview, citation_search, clinvar_annotation_health_check
        from annotation.signals.manual_signals import gene_counts_classification_withdraw_handler, \
            gene_counts_classification_publish_handler
        from classification.models import Classification, classification_withdraw_signal, \
            classification_post_publish_signal
        # pylint: enable=import-outside-toplevel,unused-import

        # Entrez wants both email and API key
        if entrez_api_key := getattr(settings, "ANNOTATION_ENTREZ_API_KEY", None):
            Entrez.api_key = entrez_api_key
        if entrez_email := getattr(settings, "ANNOTATION_ENTREZ_EMAIL", None):
            Entrez.email = entrez_email

        post_save.connect(clinvar_citations_post_save_handler, sender=CachedWebResource)

        classification_withdraw_signal.connect(gene_counts_classification_withdraw_handler,
                                               sender=Classification)
        classification_post_publish_signal.connect(gene_counts_classification_publish_handler,
                                                   sender=Classification)
