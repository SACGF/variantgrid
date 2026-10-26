from django.apps import AppConfig


class VariantopediaConfig(AppConfig):
    name = 'variantopedia'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        from snpdb.search import search_registry

        # Registers receivers on import - noqa: F401 keeps the unused-import autofix from
        # silently unregistering them
        from variantopedia.signals import integration_health_check  # noqa: F401
        from variantopedia.signals.hgvs_search import search_hgvs
        from variantopedia.signals.variant_search import (
            allele_search,
            search_allele_id,
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
        # pylint: enable=import-outside-toplevel,unused-import

        search_registry.register(
            variant_cosmic_search, search_variant_locus_no_ref, search_variant_locus_with_ref, allele_search,
            variant_search_vcf, search_variant_gnomad, search_variant_variant, search_variant_symbolic,
            search_variant_db_snp, search_hgvs, search_variant_id, search_allele_id, search_variant_gene_fusion,
            search_variant_gene_copy_number, search_variant_splice_event,
        )
