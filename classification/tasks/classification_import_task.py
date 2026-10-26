import logging

import celery

from classification.classification_import import process_classification_import
from classification.models import ClassificationImport
from classification.models.classification_import_run import ClassificationImportRun
from library.log_utils import log_traceback
from snpdb.enums import ImportSource


@celery.shared_task
def process_classification_import_task(classification_import_id: int, import_source: ImportSource):
    try:
        classification_import = ClassificationImport.objects.get(pk=classification_import_id)
        process_classification_import(classification_import, import_source)
    except:
        logging.error("process_classification_import_task FAILED")
        log_traceback()


@celery.shared_task
def close_abandoned_classification_import_runs_task():
    """ An import whose client stopped sending leaves its run ONGOING, which defers grouping, overlap and
        discordance work until the next import comes along to close it - close it here instead """
    ClassificationImportRun.cleanup()
