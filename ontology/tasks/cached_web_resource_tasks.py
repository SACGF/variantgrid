from annotation.tasks.cached_web_resource_tasks import CachedWebResourceTask
from ontology.gencc import store_gencc_from_web
from variantgrid.celery import app


class ClinGenCCWebResourceTask(CachedWebResourceTask):

    def _load_cached_web_resource(self, cached_web_resource):
        store_gencc_from_web(cached_web_resource)


ClinGenCCWebResourceTask = app.register_task(ClinGenCCWebResourceTask())
