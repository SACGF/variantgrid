from annotation.tasks.cached_web_resource_tasks import CachedWebResourceTask
from genes.gnomad_gene_constraint import store_gnomad_gene_constraint_from_web
from genes.models import PanelAppServer
from genes.panel_app import store_panel_app_panels_from_web
from genes.pfam import store_pfam_from_web
from variantgrid.celery import app


class PanelAppEnglandPanelsWebResourceTask(CachedWebResourceTask):
    def _load_cached_web_resource(self, cached_web_resource):
        panel_app_server = PanelAppServer.england_instance()
        store_panel_app_panels_from_web(panel_app_server, cached_web_resource)


class PanelAppAustraliaPanelsWebResourceTask(CachedWebResourceTask):
    def _load_cached_web_resource(self, cached_web_resource):
        panel_app_server = PanelAppServer.australia_instance()
        store_panel_app_panels_from_web(panel_app_server, cached_web_resource)


class GnomADGeneConstraintWebResourceTask(CachedWebResourceTask):
    def _load_cached_web_resource(self, cached_web_resource):
        store_gnomad_gene_constraint_from_web(cached_web_resource)


class PfamWebResourceTask(CachedWebResourceTask):
    def _load_cached_web_resource(self, cached_web_resource):
        store_pfam_from_web(cached_web_resource)


PanelAppEnglandPanelsWebResourceTask = app.register_task(PanelAppEnglandPanelsWebResourceTask())  # @UndefinedVariable
PanelAppAustraliaPanelsWebResourceTask = app.register_task(PanelAppAustraliaPanelsWebResourceTask())  # @UndefinedVariable
GnomADGeneConstraintWebResourceTask = app.register_task(GnomADGeneConstraintWebResourceTask())  # @UndefinedVariable
PfamWebResourceTask = app.register_task(PfamWebResourceTask())  # @UndefinedVariable