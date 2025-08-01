import logging

from django.conf import settings

from annotation.models.models import CachedWebResource
from genes.models import GeneListCategory
from genes.tasks.cached_web_resource_tasks import PanelAppEnglandPanelsWebResourceTask, \
    PanelAppAustraliaPanelsWebResourceTask, GnomADGeneConstraintWebResourceTask, PfamWebResourceTask, \
    UniProtWebResourceTask, RefSeqGeneSummaryWebResourceTask, HGNCWebResourceTask, LRGRefSeqGeneWebResourceTask, \
    RefSeqGeneInfoWebResourceTask, RefSeqSequenceInfoWebResourceTask, MANEWebResourceTask, \
    RefSeqGenePubMedCountWebResourceTask

# Remember to connect these handlers up in apps.GenesConfig
# For some reason this doesn't work as a variable, has to be stored here...
gnomad_gene_constraint_post_save_handler = CachedWebResource.named_handler_factory(settings.CACHED_WEB_RESOURCE_GNOMAD_GENE_CONSTRAINT,
                                                                                   GnomADGeneConstraintWebResourceTask)

hgnc_post_save_handler = CachedWebResource.named_handler_factory(settings.CACHED_WEB_RESOURCE_HGNC, HGNCWebResourceTask)

lrg_ref_seq_gene_post_save_handler = CachedWebResource.named_handler_factory(settings.CACHED_WEB_RESOURCE_LRG_REF_SEQ_GENE,
                                                                             LRGRefSeqGeneWebResourceTask)

mane_post_save_handler = CachedWebResource.named_handler_factory(settings.CACHED_WEB_RESOURCE_MANE,
                                                                 MANEWebResourceTask)

panel_app_england_panels_post_save_handler = CachedWebResource.named_handler_factory(settings.CACHED_WEB_RESOURCE_PANEL_APP_ENGLAND_PANELS,
                                                                                     PanelAppEnglandPanelsWebResourceTask)

panel_app_australia_panels_post_save_handler = CachedWebResource.named_handler_factory(settings.CACHED_WEB_RESOURCE_PANEL_APP_AUSTRALIA_PANELS,
                                                                                       PanelAppAustraliaPanelsWebResourceTask)

pfam_post_save_handler = CachedWebResource.named_handler_factory(settings.CACHED_WEB_RESOURCE_PFAM, PfamWebResourceTask)

refseq_gene_summary_post_save_handler = CachedWebResource.named_handler_factory(settings.CACHED_WEB_RESOURCE_REFSEQ_GENE_SUMMARY, RefSeqGeneSummaryWebResourceTask)

refseq_gene_info_post_save_handler = CachedWebResource.named_handler_factory(settings.CACHED_WEB_RESOURCE_REFSEQ_GENE_INFO, RefSeqGeneInfoWebResourceTask)

refseq_sequence_info_post_save_handler = CachedWebResource.named_handler_factory(settings.CACHED_WEB_RESOURCE_REFSEQ_SEQUENCE_INFO, RefSeqSequenceInfoWebResourceTask)

refseq_gene_pub_med_count_post_save_handler = CachedWebResource.named_handler_factory(settings.CACHED_WEB_RESOURCE_REFSEQ_GENE_PUBMED_COUNTS, RefSeqGenePubMedCountWebResourceTask)

uniprot_post_save_handler = CachedWebResource.named_handler_factory(settings.CACHED_WEB_RESOURCE_UNIPROT, UniProtWebResourceTask)

# Remember to connect these handlers up in apps.GenesConfig


def cached_third_part_gene_list_pre_delete_handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
    logging.info("Deleting all associated gene lists.")
    GeneListCategory.objects.filter(company=instance.company).delete()
