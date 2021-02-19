from django.apps import AppConfig
from django.db.models.signals import post_save, pre_delete


class GenesConfig(AppConfig):
    name = 'genes'

    def ready(self):
        from annotation.models.models import CachedWebResource
        from genes.models import CachedThirdPartyGeneList
        from genes.signals import hgnc_post_save_handler, panel_app_england_panels_post_save_handler, \
            panel_app_australia_panels_post_save_handler, pfam_post_save_handler, \
            gnomad_gene_constraint_post_save_handler, cached_third_part_gene_list_pre_delete_handler, \
            refseq_gene_summary_post_save_handler, uniprot_post_save_handler

        post_save.connect(gnomad_gene_constraint_post_save_handler, sender=CachedWebResource)
        post_save.connect(hgnc_post_save_handler, sender=CachedWebResource)
        post_save.connect(panel_app_england_panels_post_save_handler, sender=CachedWebResource)
        post_save.connect(panel_app_australia_panels_post_save_handler, sender=CachedWebResource)
        post_save.connect(pfam_post_save_handler, sender=CachedWebResource)
        post_save.connect(refseq_gene_summary_post_save_handler, sender=CachedWebResource)
        post_save.connect(uniprot_post_save_handler, sender=CachedWebResource)

        pre_delete.connect(cached_third_part_gene_list_pre_delete_handler, CachedThirdPartyGeneList)
