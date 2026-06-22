from django.apps import AppConfig
from django.db.models.signals import post_delete, post_save, pre_delete


class AnalysisConfig(AppConfig):
    name = 'analysis'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        # imported to activate receivers
        from analysis.signals import analysis_search, analysis_health_check

        from analysis.models import VariantTag
        from analysis.signals.signal_handlers import variant_tag_create, variant_tag_delete, \
            handle_vcf_import_success, handle_active_sample_gene_list_created
        from analysis.signals.source_data_invalidation import handle_sample_pre_delete, \
            handle_cohort_pre_delete, handle_trio_pre_delete, handle_pedigree_pre_delete, \
            handle_quad_pre_delete
        from genes.models import ActiveSampleGeneList
        from pedigree.models import Pedigree
        from snpdb.models import Cohort, Quad, Sample, Trio
        from upload.signals.signals import vcf_import_success_signal
        # pylint: enable=import-outside-toplevel,unused-import

        post_save.connect(variant_tag_create, sender=VariantTag)
        post_delete.connect(variant_tag_delete, sender=VariantTag)
        vcf_import_success_signal.connect(handle_vcf_import_success)
        post_save.connect(handle_active_sample_gene_list_created, sender=ActiveSampleGeneList)

        # Bump analysis source-node versions when their input data is deleted, so cached
        # q-dicts referencing now-missing CohortGenotype annotations are invalidated.
        # See SACGF/variantgrid_com#22.
        pre_delete.connect(handle_sample_pre_delete, sender=Sample)
        pre_delete.connect(handle_cohort_pre_delete, sender=Cohort)
        pre_delete.connect(handle_trio_pre_delete, sender=Trio)
        pre_delete.connect(handle_pedigree_pre_delete, sender=Pedigree)
        pre_delete.connect(handle_quad_pre_delete, sender=Quad)
