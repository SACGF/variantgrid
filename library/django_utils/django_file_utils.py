"""
Where an import pipeline's scratch files go: get_import_processing_dir / get_import_processing_filename
under settings.PRIVATE_DATA_ROOT, keyed on the pipeline pk.

get_import_processing_dir creates the directory, so it is for writers only. Anything that just wants to
name a directory - to remove it, or to test whether a path is under it - uses import_processing_dir_path,
otherwise merely asking where a thing would be leaves an empty directory behind (#928).
"""
import logging
import os
import shutil

from django.conf import settings

from library.utils.file_utils import mk_path


def import_processing_dir_path(pk, prefix='pipeline') -> str:
    """ Where this pk's scratch dir is (or would be) - creates nothing """
    return os.path.join(settings.IMPORT_PROCESSING_DIR, f"{prefix}_{pk}")


def get_import_processing_dir(pk, prefix='pipeline') -> str:
    upd = import_processing_dir_path(pk, prefix)
    mk_path(upd)
    return upd


def get_import_processing_filename(pk, base_filename, prefix='pipeline') -> str:
    processing_dir = get_import_processing_dir(pk, prefix)
    filename = os.path.join(processing_dir, base_filename)
    return filename


def remove_import_processing_dir(pk, prefix='pipeline'):
    """ Best effort - a dir that's already gone (eg a retry cleaned it up) is not an error """
    import_processing_dir = import_processing_dir_path(pk, prefix)
    if os.path.exists(import_processing_dir):
        logging.info("Removing import processing dir: '%s'", import_processing_dir)
        shutil.rmtree(import_processing_dir, ignore_errors=True)
