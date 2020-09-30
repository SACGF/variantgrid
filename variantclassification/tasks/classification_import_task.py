import logging

import celery

from library.log_utils import log_traceback
from variantclassification.models import VariantClassificationImport
from variantclassification.variant_classification_import import process_variant_classification_import


@celery.task
def process_variant_classification_import_task(variant_classification_import_id, import_source):
    try:
        variant_classification_import = VariantClassificationImport.objects.get(pk=variant_classification_import_id)
        process_variant_classification_import(variant_classification_import, import_source)
    except:
        logging.error("process_variant_classification_import_task FAILED")
        log_traceback()