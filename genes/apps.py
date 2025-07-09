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
            refseq_sequence_info_post_save_handler, refseq_gene_pub_med_count_post_save_handler, \
            uniprot_post_save_handler
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
        post_save.connect(refseq_gene_pub_med_count_post_save_handler, sender=CachedWebResource)
        post_save.connect(uniprot_post_save_handler, sender=CachedWebResource)

        pre_delete.connect(cached_third_part_gene_list_pre_delete_handler, CachedThirdPartyGeneList)
