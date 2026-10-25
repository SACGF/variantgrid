"""
Custom model fields. `UTCTimestampField` is a DateTimeField stored in a PostgreSQL `timestamp` (without time zone)
column holding UTC - the ClassificationModification created/modified columns, kept that type so existing tables need no
rewrite.
"""
import datetime
from zoneinfo import ZoneInfo

from django.conf import settings
from django.db import models
from django.utils import timezone

UTC = ZoneInfo("UTC")


class UTCTimestampField(models.DateTimeField):
    """ A DateTimeField in a time zone naive `timestamp` column: written as UTC, read back as aware UTC """

    description = "Date (with time and no time zone)"

    def db_type(self, connection):
        return "timestamp"

    def from_db_value(self, value, expression, connection):  # pylint: disable=unused-argument
        if isinstance(value, datetime.datetime):
            if settings.USE_TZ and timezone.is_naive(value):
                value = value.replace(tzinfo=UTC)
            return value
        return self.to_python(value)

    def get_prep_value(self, value):
        if isinstance(value, datetime.datetime):
            value = value.astimezone(UTC)
        return super().get_prep_value(value)
