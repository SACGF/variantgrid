from django.apps import AppConfig
from django.db.models.signals import post_delete, post_save

class AnalysisConfig(AppConfig):
    name = 'analysis'

    def ready(self):
        from analysis.models import VariantTag
        from analysis.signals.signal_handlers import variant_tag_create, variant_tag_delete

        post_save.connect(variant_tag_create, sender=VariantTag)
        post_delete.connect(variant_tag_delete, sender=VariantTag)
