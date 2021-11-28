import logging

import celery

from classification.classification_import import process_classification_import
from classification.models import ClassificationImport
from library.log_utils import log_traceback


@celery.shared_task
def process_classification_import_task(classification_import_id, import_source):
    try:
        classification_import = ClassificationImport.objects.get(pk=classification_import_id)
        process_classification_import(classification_import, import_source)
    except:
        logging.error("process_classification_import_task FAILED")
        log_traceback()