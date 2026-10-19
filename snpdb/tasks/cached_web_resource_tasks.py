"""
CachedWebResourceTask: the celery task base for loading a CachedWebResource. It marks the resource
IMPORTING, calls _load_cached_web_resource, then records SUCCESS or ERROR (with the traceback in description).
Subclasses live with the data they load (annotation, genes and ontology tasks/cached_web_resource_tasks.py).
"""
import logging

from celery.app.task import Task

from library.log_utils import get_traceback
from snpdb.models.models_cached_web_resource import CachedWebResource
from snpdb.models.models_enums import ImportStatus


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
        except Exception:
            tb = get_traceback()
            cached_web_resource.import_status = ImportStatus.ERROR
            cached_web_resource.description = tb
            logging.error(tb)

        cached_web_resource.save()
