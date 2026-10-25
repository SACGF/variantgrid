""" Who may write sequencing data - the pipeline's API user. Used by the seqauto REST write views
    (seqauto/views_rest.py:SeqAutoWritePermission) and by uploads that write onto seqauto records
    (upload/models/models.py:FileUpload.check_can_write_sequencing_data). """

from django.conf import settings
from django.contrib.auth.models import User


def can_write_seqauto(user: User) -> bool:
    return user.is_superuser or user.groups.filter(name=settings.SEQAUTO_API_WRITE_GROUP).exists()
