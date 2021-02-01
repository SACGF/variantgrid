from annotation.tasks.cached_web_resource_tasks import CachedWebResourceTask
from genes.gnomad_gene_constraint import store_gnomad_gene_constraint_from_web
from genes.hgnc import store_hgnc_from_web
from genes.models import PanelAppServer
from genes.panel_app import store_panel_app_panels_from_web
from genes.pfam import store_pfam_from_web
from genes.refseq import store_refseq_gene_summary_from_web
from genes.rvis import store_rvis_from_web
from genes.uniprot import store_uniprot_from_web
from variantgrid.celery import app


class GnomADGeneConstraintWebResourceTask(CachedWebResourceTask):
    def _load_cached_web_resource(self, cached_web_resource):
        store_gnomad_gene_constraint_from_web(cached_web_resource)


class HGNCWebResourceTask(CachedWebResourceTask):
    def _load_cached_web_resource(self, cached_web_resource):
        store_hgnc_from_web(cached_web_resource)


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


class RVISWebResourceTask(CachedWebResourceTask):
    def _load_cached_web_resource(self, cached_web_resource):
        store_rvis_from_web(cached_web_resource)


class UniProtWebResourceTask(CachedWebResourceTask):
    def _load_cached_web_resource(self, cached_web_resource):
        store_uniprot_from_web(cached_web_resource)


GnomADGeneConstraintWebResourceTask = app.register_task(GnomADGeneConstraintWebResourceTask())
HGNCWebResourceTask = app.register_task(HGNCWebResourceTask())
PanelAppEnglandPanelsWebResourceTask = app.register_task(PanelAppEnglandPanelsWebResourceTask())
PanelAppAustraliaPanelsWebResourceTask = app.register_task(PanelAppAustraliaPanelsWebResourceTask())
PfamWebResourceTask = app.register_task(PfamWebResourceTask())
RefSeqGeneSummaryWebResourceTask = app.register_task(RefSeqGeneSummaryWebResourceTask())
RVISWebResourceTask = app.register_task(RVISWebResourceTask())
UniProtWebResourceTask = app.register_task(UniProtWebResourceTask())
