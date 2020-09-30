from celery.app.task import Task
import logging

from annotation.clingen.gene_validity import store_clingen_gene_validity_curations_from_web
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


class ClinGenValidityCurationsWebResourceTask(CachedWebResourceTask):

    def _load_cached_web_resource(self, cached_web_resource):
        store_clingen_gene_validity_curations_from_web(cached_web_resource)


ClinGenValidityCurationsWebResourceTask = app.register_task(ClinGenValidityCurationsWebResourceTask())  # @UndefinedVariable
