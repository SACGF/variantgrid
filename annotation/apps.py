from django.apps.config import AppConfig
from django.db import ProgrammingError
from django.db.models.signals import post_save


class AnnotationConfig(AppConfig):
    name = 'annotation'

    def ready(self):
        from django.conf import settings
        from annotation.models import CachedWebResource
        from annotation.signals import clingen_post_save_handler

        # Entrez wants both email and API key
        if entrez_api_key := getattr(settings, "ANNOTATION_ENTREZ_API_KEY", None):
            from Bio import Entrez
            Entrez.api_key = entrez_api_key
        if entrez_email := getattr(settings, "ANNOTATION_ENTREZ_EMAIL", None):
            from Bio import Entrez
            Entrez.email = entrez_email

        post_save.connect(clingen_post_save_handler, sender=CachedWebResource)
        try:
            GeneCountType = self.get_model('GeneCountType')
            if GeneCountType.objects.filter(enabled=True, uses_classifications=True).exists():
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
