from django.apps import AppConfig
from django.db.models.signals import post_delete, post_save, pre_delete


class AnalysisConfig(AppConfig):
    name = 'analysis'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        # imported to activate receivers

        from analysis.models import Analysis, VariantTag

        # Registers receivers on import - noqa: F401 keeps the unused-import autofix from
        # silently unregistering them
        from analysis.signals import analysis_health_check, analysis_search  # noqa: F401
        from analysis import user_awards  # noqa: F401  # registers award definitions on import
        from analysis.signals.signal_handlers import (
            analysis_pre_delete,
            handle_active_sample_gene_list_created,
            handle_vcf_import_success,
            variant_tag_create,
            variant_tag_delete,
        )
        from analysis.signals.source_data_invalidation import (
            handle_cohort_pre_delete,
            handle_duo_pre_delete,
            handle_pedigree_pre_delete,
            handle_quad_pre_delete,
            handle_sample_pre_delete,
            handle_trio_pre_delete,
        )
        from genes.models import ActiveSampleGeneList
        from pedigree.models import Pedigree
        from snpdb.models import Cohort, Duo, Quad, Sample, Trio
        from upload.signals.signals import vcf_import_success_signal
        # pylint: enable=import-outside-toplevel,unused-import

        post_save.connect(variant_tag_create, sender=VariantTag)
        post_delete.connect(variant_tag_delete, sender=VariantTag)
        vcf_import_success_signal.connect(handle_vcf_import_success)
        post_save.connect(handle_active_sample_gene_list_created, sender=ActiveSampleGeneList)

        # Stop any node still loading before the delete cascades into its caches - see
        # analysis_pre_delete. SACGF/variantgrid_com#2
        pre_delete.connect(analysis_pre_delete, sender=Analysis)

        # Bump analysis source-node versions when their input data is deleted, so cached
        # q-dicts referencing now-missing CohortGenotype annotations are invalidated.
        # See SACGF/variantgrid_com#22.
        pre_delete.connect(handle_sample_pre_delete, sender=Sample)
        pre_delete.connect(handle_cohort_pre_delete, sender=Cohort)
        pre_delete.connect(handle_trio_pre_delete, sender=Trio)
        pre_delete.connect(handle_pedigree_pre_delete, sender=Pedigree)
        pre_delete.connect(handle_quad_pre_delete, sender=Quad)
        pre_delete.connect(handle_duo_pre_delete, sender=Duo)
