from django.apps import AppConfig
from django.db.models.signals import post_delete, post_save


class AnalysisConfig(AppConfig):
    name = 'analysis'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel
        from analysis.models import VariantTag
        from analysis.signals.signal_handlers import variant_tag_create, variant_tag_delete

        # imported to activate receivers
        from analysis.signals import analysis_health_check  # pylint: disable=unused-import
        from analysis.signals import analysis_search
        # pylint: enable=import-outside-toplevel

        post_save.connect(variant_tag_create, sender=VariantTag)
        post_delete.connect(variant_tag_delete, sender=VariantTag)
