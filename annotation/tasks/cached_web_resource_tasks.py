import logging

from celery.app.task import Task

from annotation.clingen.clinvar_citations import store_clinvar_citations_from_web
from annotation.models.models import CachedWebResource
from library.log_utils import get_traceback
from snpdb.models.models_enums import ImportStatus
from variantgrid.celery import app


class CachedWebResourceTask(Task):
    abstract = True

    def _load_cached_web_resource(self, cached_web_resource):
        pass

    def run(self, cached_web_resource_id):
        cached_web_resource = CachedWebResource.objects.get(pk=cached_web_resource_id)
        cached_web_resource.import_status = ImportStatus.IMPORTING
        cached_web_resource.save()
        try:

            self._load_cached_web_resource(cached_web_resource)
            cached_web_resource.import_status = ImportStatus.SUCCESS
        except:
            tb = get_traceback()
            cached_web_resource.import_status = ImportStatus.ERROR
            cached_web_resource.description = tb
            logging.error(tb)

        cached_web_resource.save()


class ClinVarCitationsWebResourceTask(CachedWebResourceTask):

    def _load_cached_web_resource(self, cached_web_resource):
        store_clinvar_citations_from_web(cached_web_resource)


ClinVarCitationsWebResourceTask = app.register_task(ClinVarCitationsWebResourceTask())
