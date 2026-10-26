from annotation.clingen.clinvar_citations import store_clinvar_citations_from_web
from snpdb.tasks.cached_web_resource_tasks import CachedWebResourceTask
from variantgrid.celery import app


class ClinVarCitationsWebResourceTask(CachedWebResourceTask):

    def _load_cached_web_resource(self, cached_web_resource):
        store_clinvar_citations_from_web(cached_web_resource)


ClinVarCitationsWebResourceTask = app.register_task(ClinVarCitationsWebResourceTask())
