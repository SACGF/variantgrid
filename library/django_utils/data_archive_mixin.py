from django.contrib.auth.models import User
from django.db import models
from django.db.models.deletion import SET_NULL


class DataArchivedError(Exception):
    """ Raised when code tries to read underlying data that has been archived.

        Callers either propagate to a UI surface (banner / badge) or catch in
        a source-node helper that maps to NodeStatus.ERROR_CONFIGURATION.
    """

    def __init__(self, obj):
        self.obj = obj
        date = getattr(obj, "data_archived_date", None)
        by = getattr(obj, "data_archived_by", None)
        reason = getattr(obj, "data_archive_reason", None) or ""
        date_str = date.strftime("%Y-%m-%d") if date else "unknown date"
        by_str = str(by) if by else "unknown user"
        msg = f"Source data archived on {date_str} by {by_str}"
        if reason:
            msg += f": {reason}"
        super().__init__(msg)


class DataArchiveMixin(models.Model):
    """ Stamps a model when its underlying data has been dropped but the row is kept.

        @see https://github.com/SACGF/variantgrid/issues/1536
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
