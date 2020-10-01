import logging

import celery

from library.log_utils import log_traceback
from classification.models import ClassificationImport
from classification.classification_import import process_classification_import


@celery.task
def process_classification_import_task(classification_import_id, import_source):
    try:
        classification_import = ClassificationImport.objects.get(pk=classification_import_id)
        process_classification_import(classification_import, import_source)
    except:
        logging.error("process_classification_import_task FAILED")
        log_traceback()