import logging

from importlib import metadata
from django.apps import AppConfig
from django.db.models.signals import post_save, pre_delete


class GenesConfig(AppConfig):
    name = 'genes'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        from genes.signals import gene_search, gene_symbol_search, transcript_search

        from annotation.models.models import CachedWebResource
        from genes.models import CachedThirdPartyGeneList
        from genes.signals.manual_signals import hgnc_post_save_handler, lrg_ref_seq_gene_post_save_handler, \
            mane_post_save_handler, \
            panel_app_england_panels_post_save_handler, panel_app_australia_panels_post_save_handler, \
            pfam_post_save_handler, \
            gnomad_gene_constraint_post_save_handler, cached_third_part_gene_list_pre_delete_handler, \
            refseq_gene_summary_post_save_handler, refseq_gene_info_post_save_handler, \
            refseq_sequence_info_post_save_handler, uniprot_post_save_handler
        # pylint: enable=import-outside-toplevel,unused-import

        post_save.connect(gnomad_gene_constraint_post_save_handler, sender=CachedWebResource)
        post_save.connect(hgnc_post_save_handler, sender=CachedWebResource)
        post_save.connect(lrg_ref_seq_gene_post_save_handler, sender=CachedWebResource)
        post_save.connect(mane_post_save_handler, sender=CachedWebResource)
        post_save.connect(panel_app_england_panels_post_save_handler, sender=CachedWebResource)
        post_save.connect(panel_app_australia_panels_post_save_handler, sender=CachedWebResource)
        post_save.connect(pfam_post_save_handler, sender=CachedWebResource)
        post_save.connect(refseq_gene_summary_post_save_handler, sender=CachedWebResource)
        post_save.connect(refseq_gene_info_post_save_handler, sender=CachedWebResource)
        post_save.connect(refseq_sequence_info_post_save_handler, sender=CachedWebResource)
        post_save.connect(uniprot_post_save_handler, sender=CachedWebResource)

        pre_delete.connect(cached_third_part_gene_list_pre_delete_handler, CachedThirdPartyGeneList)

        # Check library versions to make sure bug fixes have been applied

        def _test_biocommons_hgvs():
            from snpdb.models import GenomeBuild
            from genes.hgvs import HGVSMatcher

            try:
                matcher = HGVSMatcher(GenomeBuild.grch38())
                matcher.get_variant_coordinate("NC_000006.12:g.49949407_49949408=")
            except TypeError:
                logging.error("Your Biocommons HGVS library is out of date, please run: " +
                              "'python3 -m pip install git+https://github.com/davmlaw/hgvs@variantgrid#egg=hgvs'")

        MINIMUM_VERSIONS = {
            "cdot": (0, 2, 21),
            "hgvs": _test_biocommons_hgvs,
            "pyhgvs": (0, 12, 4),
        }

        for name, version_required in MINIMUM_VERSIONS.items():
            if callable(version_required):
                version_required()
            else:
                version_str = metadata.version(name)
                version = tuple([int(i) for i in version_str.split(".")])
                if version < version_required:
                    logging.error("Library %s (%s) requires version >= %s", name, version, version_required)
