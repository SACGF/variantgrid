from django.apps import AppConfig
from django.db.models.signals import post_delete, post_save


class AnalysisConfig(AppConfig):
    name = 'analysis'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        # imported to activate receivers
        from analysis.signals import analysis_search, analysis_health_check

        from analysis.models import VariantTag
        from analysis.signals.signal_handlers import variant_tag_create, variant_tag_delete, vcf_import_success
        from snpdb.models import VCF
        from upload.signals.signals import vcf_import_success_signal
        # pylint: enable=import-outside-toplevel,unused-import

        post_save.connect(variant_tag_create, sender=VariantTag)
        post_delete.connect(variant_tag_delete, sender=VariantTag)
        vcf_import_success_signal.connect(vcf_import_success)
