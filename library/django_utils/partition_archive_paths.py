import os

from django.conf import settings


def get_partition_archive_path(filename: str) -> str:
    """ Resolve the on-disk path for a partition dump file.

        Files live under <PARTITION_ARCHIVE_DIR>/<postgres_db_name>/<filename>
        so test/prod/demo deploys on shared storage do not clobber each other.
    """
    db_name = settings.DATABASES["default"]["NAME"]
    base = os.path.join(settings.PARTITION_ARCHIVE_DIR, db_name)
    os.makedirs(base, exist_ok=True)
    return os.path.join(base, filename)
