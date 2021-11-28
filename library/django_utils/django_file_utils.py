import os

from django.conf import settings

from library.file_utils import mk_path


def get_import_processing_dir(pk, prefix='pipeline'):
    output_dir = f"{prefix}_{pk}"
    upd = os.path.join(settings.IMPORT_PROCESSING_DIR, output_dir)
    mk_path(upd)
    return upd


def get_import_processing_filename(pk, base_filename, prefix='pipeline'):
    processing_dir = get_import_processing_dir(pk, prefix)
    filename = os.path.join(processing_dir, base_filename)
    return filename
