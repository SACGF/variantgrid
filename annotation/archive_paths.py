import os

from django.conf import settings


def get_annotation_archive_path(partition_filename: str) -> str:
    """ Resolve the on-disk path for an annotation-partition dump file.

        The DB name is included in the namespace so test/prod/demo deploys on
        shared NFS/EFS don't clobber each other.

        @see claude/issue_1536_data_archive_plan.md (consumed by #1537)
    """
    db_name = settings.DATABASES["default"]["NAME"]
    base = os.path.join(settings.ANNOTATION_ARCHIVE_DIR, db_name)
    os.makedirs(base, exist_ok=True)
    return os.path.join(base, partition_filename)
