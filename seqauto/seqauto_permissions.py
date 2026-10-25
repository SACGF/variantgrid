""" Who may write sequencing data - the pipeline's API user. Used by the seqauto REST write views
    (seqauto/views_rest.py:SeqAutoWritePermission) and when an upload is linked to sequencing data
    (upload/models/models.py:FileUpload.check_may_link_to_sequencing_data). """

from django.conf import settings
from django.contrib.auth.models import User

from library.guardian_utils import user_in_group_or_superuser


def can_write_seqauto(user: User) -> bool:
    return user_in_group_or_superuser(user, settings.SEQAUTO_API_WRITE_GROUP)
