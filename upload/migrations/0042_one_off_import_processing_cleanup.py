import os

from django.conf import settings
from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _has_import_processing_files(apps):  # pylint: disable=unused-argument
    """ #928: several paths never removed their import_processing scratch, so an existing deployment has
        years of it. Only worth surfacing where there is something to reclaim """
    import_processing_dir = settings.IMPORT_PROCESSING_DIR
    try:
        return bool(os.listdir(import_processing_dir))
    except OSError:
        return False


class Migration(migrations.Migration):
    dependencies = [
        ("upload", "0041_alter_fileupload_file_type_and_more"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["import_processing_cleanup"]),
                        note="Reclaim import_processing scratch left behind by imports that never "
                             "cleaned up (#928). Run with --dry-run first to see what it would remove",
                        test=_has_import_processing_files),
    ]
