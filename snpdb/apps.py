import logging
import mimetypes
import sys

from django.apps import AppConfig
from django.conf import settings
from django.db.models.signals import post_save


class SnpdbConfig(AppConfig):
    name = 'snpdb'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        from django.contrib.auth.models import Group, User

        from snpdb import checks  # noqa: F401  # registers system checks on import
        from snpdb import user_award_definitions  # noqa: F401  # registers award definitions on import
        from snpdb.models import Trio

        # Registers receivers on import - noqa: F401 keeps the unused-import autofix from
        # silently unregistering them
        from snpdb.signals import (  # noqa: F401
            common_variants_classification_changed,
            disk_usage_health_check,
            jobs_autopause,  # registers worker_ready crash-safety auto-pause
            somalier_patient_relate,
            variant_zygosity_preview_extra,
            vcf_health_check,
        )
        from snpdb.search import search_registry
        from snpdb.signals.clinvar_export_search import clinvar_export_batch_search, clinvar_id_search
        from snpdb.signals.cohort_search import search_cohort
        from snpdb.signals.duo_search import search_duo
        from snpdb.signals.genomics_search import contig_search, genome_build_search
        from snpdb.signals.lab_search import lab_search
        from snpdb.signals.organization_search import organization_search
        from snpdb.signals.quad_search import search_quad
        from snpdb.signals.sample_search import sample_search
        from snpdb.signals.scv_search import scv_search
        from snpdb.signals.signal_handlers import (
            group_post_save_handler,
            trio_post_save_handler,
            user_post_save_handler,
        )
        from snpdb.signals.trio_search import search_trio
        from snpdb.signals.user_search import user_search
        from snpdb.signals.variant_search import (
            allele_search,
            search_allele_id,
            search_hgvs,
            search_variant_db_snp,
            search_variant_gene_copy_number,
            search_variant_gene_fusion,
            search_variant_gnomad,
            search_variant_id,
            search_variant_locus_no_ref,
            search_variant_locus_with_ref,
            search_variant_splice_event,
            search_variant_symbolic,
            search_variant_variant,
            variant_cosmic_search,
            variant_search_vcf,
        )
        from snpdb.signals.vcf_search import vcf_search
        # pylint: enable=import-outside-toplevel,unused-import

        search_registry.register(
            clinvar_id_search, clinvar_export_batch_search,
            search_cohort,
            search_duo,
            genome_build_search, contig_search,
            lab_search,
            organization_search,
            search_quad,
            sample_search,
            scv_search,
            search_trio,
            user_search,
            variant_cosmic_search, search_variant_locus_no_ref, search_variant_locus_with_ref, allele_search,
            variant_search_vcf, search_variant_gnomad, search_variant_variant, search_variant_symbolic,
            search_variant_db_snp, search_hgvs, search_variant_id, search_allele_id, search_variant_gene_fusion,
            search_variant_gene_copy_number, search_variant_splice_event,
            vcf_search,
        )

        if not settings.UNIT_TEST:
            # Add newly created users to public group
            post_save.connect(user_post_save_handler, sender=User)

        if settings.REQUESTS_DISABLE_IPV6:
            import requests  # pylint: disable=import-outside-toplevel
            requests.packages.urllib3.util.connection.HAS_IPV6 = False

        # Make global settings share read only with this group by default
        post_save.connect(group_post_save_handler, sender=Group)

        # Add newly created users to public group
        post_save.connect(trio_post_save_handler, sender=Trio)

        # Disable annoying matplotlib findfont messages
        logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

        if sys.version_info < (3, 10):
            raise SystemExit("VariantGrid requires Python 3.10 or later.")

        # So static serve of MEDIA_ROOT just prompts to download VCF (used to rename to vcf.vcf)
        mimetypes.add_type("application/octet-stream", ".vcf", strict=True)
