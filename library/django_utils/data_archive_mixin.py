from django.contrib.auth.models import User
from django.db import models
from django.db.models.deletion import SET_NULL


class DataArchiveMixin(models.Model):
    """ Stamps a model when its underlying data has been dropped but the row is kept.

        @see claude/issue_1536_data_archive_plan.md
    """
    data_archived_date = models.DateTimeField(null=True, blank=True)
    data_archived_by = models.ForeignKey(
        User, null=True, blank=True, on_delete=SET_NULL, related_name="+"
    )
    data_archive_reason = models.TextField(null=True, blank=True)
    data_restorable_from = models.TextField(null=True, blank=True)

    class Meta:
        abstract = True

    @property
    def data_archived(self) -> bool:
        return self.data_archived_date is not None
