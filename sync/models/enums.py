from django.db.models import TextChoices


class SyncStatus(TextChoices):
    IN_PROGRESS = 'P', 'In Progress'
    SUCCESS = 'S', 'Success'
    # With no records, we had nothing to upload or download
    # there was no failure, but no success
    NO_RECORDS = 'N', 'No Records'
    FAILED = 'F', 'Failed'
