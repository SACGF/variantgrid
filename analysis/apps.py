from django.apps import AppConfig
from django.db.models.signals import post_delete, post_save

class AnalysisConfig(AppConfig):
    name = 'analysis'

    def ready(self):
        from analysis.models import VariantTag
        from analysis.signals.signal_handlers import variant_tag_create, variant_tag_delete, \
            handle_vcf_import_success, handle_active_sample_gene_list_created
        from genes.models import ActiveSampleGeneList
        from upload.signals.signals import vcf_import_success_signal

        post_save.connect(variant_tag_create, sender=VariantTag)
        post_delete.connect(variant_tag_delete, sender=VariantTag)
        vcf_import_success_signal.connect(handle_vcf_import_success)
        post_save.connect(handle_active_sample_gene_list_created, sender=ActiveSampleGeneList)
