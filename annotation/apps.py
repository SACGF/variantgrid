from django.apps.config import AppConfig
from django.db.models.signals import post_delete, post_save


class AnnotationConfig(AppConfig):
    name = 'annotation'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        # imported to activate receivers

        from Bio import Entrez
        from django.conf import settings
        from annotation.models import AnnotationRun, CachedWebResource
        from annotation.signals.manual_signals import clinvar_citations_post_save_handler, \
            annotation_run_complete_signal, annotation_run_discarded_signal
        from annotation.signals.annotation_run_cleanup import annotation_run_complete_cleanup_handler, \
            annotation_run_discarded_cleanup_handler, annotation_run_post_delete_handler
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

        # #1670: every route by which an AnnotationRun's on-disk output is reclaimed
        annotation_run_complete_signal.connect(annotation_run_complete_cleanup_handler)
        annotation_run_discarded_signal.connect(annotation_run_discarded_cleanup_handler)
        post_delete.connect(annotation_run_post_delete_handler, sender=AnnotationRun)

        classification_withdraw_signal.connect(gene_counts_classification_withdraw_handler,
                                               sender=Classification)
        classification_post_publish_signal.connect(gene_counts_classification_publish_handler,
                                                   sender=Classification)
