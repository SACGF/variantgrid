from django.apps.config import AppConfig
from django.core.cache import cache
from django.db import ProgrammingError
from django.db.models.signals import post_save


def _has_classification_gene_count_type(GeneCountType):
    cache_key = "has_classification_gene_count_type"
    exists = cache.get(cache_key)
    if exists is None:
        exists = GeneCountType.objects.filter(enabled=True, uses_classifications=True).exists()
        cache.set(cache_key, exists)
    return exists



class AnnotationConfig(AppConfig):
    name = 'annotation'

    def ready(self):
        # pylint: disable=import-outside-toplevel
        from Bio import Entrez
        from django.conf import settings
        from annotation.models import CachedWebResource
        from annotation.signals import clinvar_citations_post_save_handler
        # pylint: enable=import-outside-toplevel

        # Entrez wants both email and API key
        if entrez_api_key := getattr(settings, "ANNOTATION_ENTREZ_API_KEY", None):
            Entrez.api_key = entrez_api_key
        if entrez_email := getattr(settings, "ANNOTATION_ENTREZ_EMAIL", None):
            Entrez.email = entrez_email

        post_save.connect(clinvar_citations_post_save_handler, sender=CachedWebResource)

        try:
            GeneCountType = self.get_model('GeneCountType')
            if _has_classification_gene_count_type(GeneCountType):
                from annotation.signals import gene_counts_classification_withdraw_handler, \
                    gene_counts_classification_publish_handler
                from classification.models import Classification, \
                    classification_withdraw_signal, classification_post_publish_signal

                classification_withdraw_signal.connect(gene_counts_classification_withdraw_handler,
                                                       sender=Classification)
                classification_post_publish_signal.connect(gene_counts_classification_publish_handler,
                                                           sender=Classification)
        except ProgrammingError:
            pass  # Need to allow DB migrations adding these fields to run
