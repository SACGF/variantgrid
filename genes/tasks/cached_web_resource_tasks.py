from annotation.tasks.cached_web_resource_tasks import CachedWebResourceTask
from genes.cached_web_resource.gnomad_gene_constraint import store_gnomad_gene_constraint_from_web
from genes.cached_web_resource.hgnc import store_hgnc_from_web
from genes.cached_web_resource.lrg_ref_seq_gene import store_lrg_ref_seq_gene_from_web
from genes.cached_web_resource.mane import store_mane_from_web
from genes.cached_web_resource.pfam import store_pfam_from_web
from genes.cached_web_resource.refseq import store_refseq_gene_summary_from_web, store_refseq_gene_info_from_web, \
    store_refseq_sequence_info_from_web, store_gene2pubmed_from_web
from genes.cached_web_resource.uniprot import store_uniprot_from_web
from genes.models import PanelAppServer
from genes.panel_app import store_panel_app_panels_from_web
from variantgrid.celery import app


class GnomADGeneConstraintWebResourceTask(CachedWebResourceTask):
    def _load_cached_web_resource(self, cached_web_resource):
        store_gnomad_gene_constraint_from_web(cached_web_resource)


class HGNCWebResourceTask(CachedWebResourceTask):
    def _load_cached_web_resource(self, cached_web_resource):
        store_hgnc_from_web(cached_web_resource)


class LRGRefSeqGeneWebResourceTask(CachedWebResourceTask):
    def _load_cached_web_resource(self, cached_web_resource):
        store_lrg_ref_seq_gene_from_web(cached_web_resource)


class MANEWebResourceTask(CachedWebResourceTask):
    def _load_cached_web_resource(self, cached_web_resource):
        store_mane_from_web(cached_web_resource)


class PanelAppAustraliaPanelsWebResourceTask(CachedWebResourceTask):
    def _load_cached_web_resource(self, cached_web_resource):
        panel_app_server = PanelAppServer.australia_instance()
        store_panel_app_panels_from_web(panel_app_server, cached_web_resource)


class PanelAppEnglandPanelsWebResourceTask(CachedWebResourceTask):
    def _load_cached_web_resource(self, cached_web_resource):
        panel_app_server = PanelAppServer.england_instance()
        store_panel_app_panels_from_web(panel_app_server, cached_web_resource)


class PfamWebResourceTask(CachedWebResourceTask):
    def _load_cached_web_resource(self, cached_web_resource):
        store_pfam_from_web(cached_web_resource)


class RefSeqGeneSummaryWebResourceTask(CachedWebResourceTask):
    def _load_cached_web_resource(self, cached_web_resource):
        store_refseq_gene_summary_from_web(cached_web_resource)


class RefSeqGeneInfoWebResourceTask(CachedWebResourceTask):
    def _load_cached_web_resource(self, cached_web_resource):
        store_refseq_gene_info_from_web(cached_web_resource)


class RefSeqSequenceInfoWebResourceTask(CachedWebResourceTask):
    def _load_cached_web_resource(self, cached_web_resource):
        store_refseq_sequence_info_from_web(cached_web_resource)


class RefSeqGenePubMedCountWebResourceTask(CachedWebResourceTask):
    def _load_cached_web_resource(self, cached_web_resource):
        store_gene2pubmed_from_web(cached_web_resource)


class UniProtWebResourceTask(CachedWebResourceTask):
    def _load_cached_web_resource(self, cached_web_resource):
        store_uniprot_from_web(cached_web_resource)


GnomADGeneConstraintWebResourceTask = app.register_task(GnomADGeneConstraintWebResourceTask())
HGNCWebResourceTask = app.register_task(HGNCWebResourceTask())
LRGRefSeqGeneWebResourceTask = app.register_task(LRGRefSeqGeneWebResourceTask())
MANEWebResourceTask = app.register_task(MANEWebResourceTask())
PanelAppEnglandPanelsWebResourceTask = app.register_task(PanelAppEnglandPanelsWebResourceTask())
PanelAppAustraliaPanelsWebResourceTask = app.register_task(PanelAppAustraliaPanelsWebResourceTask())
PfamWebResourceTask = app.register_task(PfamWebResourceTask())
RefSeqGeneSummaryWebResourceTask = app.register_task(RefSeqGeneSummaryWebResourceTask())
RefSeqGeneInfoWebResourceTask = app.register_task(RefSeqGeneInfoWebResourceTask())
RefSeqSequenceInfoWebResourceTask = app.register_task(RefSeqSequenceInfoWebResourceTask())
RefSeqGenePubMedCountWebResourceTask = app.register_task(RefSeqGenePubMedCountWebResourceTask())
UniProtWebResourceTask = app.register_task(UniProtWebResourceTask())
