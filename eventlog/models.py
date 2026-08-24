import inspect
import logging
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Optional, Union

from django.contrib.auth.models import Group, User
from django.contrib.auth.signals import user_logged_in
from django.db import models
from django.db.models import F
from django.db.models.deletion import SET_NULL
from django.utils import timezone
from model_utils.models import TimeStampedModel

from library.enums.log_level import LogLevel
from library.integration_status import IntegrationDirection


class ViewEvent(TimeStampedModel):
    # ViewEvent isn't an accurate name since it's also used for POST
    user = models.ForeignKey(User, null=True, on_delete=SET_NULL)
    view_name = models.TextField()
    args = models.JSONField(null=False, blank=True, default=dict)
    path = models.TextField()
    method = models.TextField()
    referer = models.TextField(null=True, blank=True)

    @property
    def is_get(self):
        return not self.method or self.method.upper() == "GET"


class Event(models.Model):
    user = models.ForeignKey(User, null=True, on_delete=SET_NULL)
    date = models.DateTimeField(null=False)
    app_name = models.TextField()
    name = models.TextField()
    details = models.TextField(null=True)
    severity = models.CharField(max_length=1, choices=LogLevel.CHOICES, default=LogLevel.INFO)
    filename = models.TextField(null=True)

    def can_write(self, user_or_group: Union[User, Group]) -> bool:
        return isinstance(user_or_group, User) and (user_or_group.is_superuser or self.user == user_or_group)

    def __str__(self):
        MAX_DETAILS_LENGTH = 200
        details = self.details
        if details and (len(details) > MAX_DETAILS_LENGTH):
            details = details[:MAX_DETAILS_LENGTH // 2] + " ... " + details[-MAX_DETAILS_LENGTH // 2:]
        user_msg = ''
        if self.user:
            user_msg = f"User: {self.user}, "

        return f"{user_msg}{self.app_name}/{self.name}: {details}"


class IntegrationActivityStatus(models.TextChoices):
    SUCCESS = "S", "Success"
    NO_CHANGE = "N", "No change"
    ERROR = "E", "Error"


@dataclass
class IntegrationActivityTracker:
    """ Handed to the body of IntegrationActivity.track() so it can say whether it wrote anything """
    changed: bool = False


class IntegrationActivity(TimeStampedModel):
    """ One row per integration, updated in place - a rolling "last seen" for external systems
        that keep no run log of their own. Bounded however chatty the integration is """
    key = models.TextField(unique=True)  # "sapath-mocha-api", "api:specimen_measure_bulk_create"
    name = models.TextField()  # "Mocha API" - the label on the Server Status row
    direction = models.CharField(max_length=1, choices=IntegrationDirection.choices,
                                 default=IntegrationDirection.INBOUND)
    last_attempt = models.DateTimeField(null=True, blank=True)
    last_success = models.DateTimeField(null=True, blank=True)
    last_change = models.DateTimeField(null=True, blank=True)  # last attempt that actually wrote something
    last_error = models.DateTimeField(null=True, blank=True)
    last_error_message = models.TextField(null=True, blank=True)
    attempt_count = models.IntegerField(default=0)
    success_count = models.IntegerField(default=0)
    error_count = models.IntegerField(default=0)

    class Meta:
        verbose_name_plural = "Integration activity"

    def __str__(self):
        return self.name or self.key

    @property
    def failing(self) -> bool:
        """ True when the most recent outcome was a failure """
        if not self.last_error:
            return False
        return not self.last_success or self.last_error > self.last_success

    @classmethod
    def _upsert(cls, key: str, name: Optional[str], direction: str, updates: dict):
        """ A single row-level UPDATE once the row exists, which is every call after the first """
        if name:
            updates["name"] = name
        updates["direction"] = direction
        # .update() bypasses auto_now, so keep modified in step with the stamps we're writing
        updates["modified"] = timezone.now()
        if not cls.objects.filter(key=key).update(**updates):
            cls.objects.get_or_create(key=key, defaults={"name": name or key, "direction": direction})
            cls.objects.filter(key=key).update(**updates)

    @classmethod
    def record(
            cls,
            key: str,
            name: Optional[str] = None,
            direction: str = IntegrationDirection.INBOUND,
            status: str = IntegrationActivityStatus.SUCCESS,
            changed: bool = False,
            message: Optional[str] = None):
        """ Records one attempt and its outcome. Counters go through F() so concurrent workers add up """
        timestamp = timezone.now()
        updates = {
            "last_attempt": timestamp,
            "attempt_count": F("attempt_count") + 1,
        }
        if status == IntegrationActivityStatus.ERROR:
            updates["last_error"] = timestamp
            updates["last_error_message"] = message
            updates["error_count"] = F("error_count") + 1
        else:
            updates["last_success"] = timestamp
            updates["success_count"] = F("success_count") + 1
            if changed:
                updates["last_change"] = timestamp
        cls._upsert(key, name, direction, updates)

    @classmethod
    @contextmanager
    def track(cls, key: str, name: Optional[str] = None, direction: str = IntegrationDirection.INBOUND):
        """
        with IntegrationActivity.track("sapath-mocha-api", name="Mocha API") as activity:
            activity.changed = bool(import_mocha_sample_extractions())
        """
        cls._upsert(key, name, direction, {
            "last_attempt": timezone.now(),
            "attempt_count": F("attempt_count") + 1,
        })
        tracker = IntegrationActivityTracker()
        try:
            yield tracker
        except Exception as e:
            timestamp = timezone.now()
            cls._upsert(key, name, direction, {
                "last_error": timestamp,
                "last_error_message": str(e),
                "error_count": F("error_count") + 1,
            })
            raise
        timestamp = timezone.now()
        updates = {
            "last_success": timestamp,
            "success_count": F("success_count") + 1,
        }
        if tracker.changed:
            updates["last_change"] = timestamp
        cls._upsert(key, name, direction, updates)


def create_login_event(sender, user, request, **kwargs):  # pylint: disable=unused-argument
    event = Event(user=user,
                  app_name='snpdb',
                  name='login',
                  date=timezone.now())
    event.save()


def create_event(
        user: Optional[User],
        name: str,
        details: Optional[str] = None,
        filename: Optional[str] = None,
        severity=LogLevel.INFO,
        app_name: Optional[str] = None,
        log=True) -> Event:

    if app_name is None:
        frm = inspect.stack()[1]
        mod = inspect.getmodule(frm[0])
        app_name = mod.__name__.split('.')[0]

    e = Event.objects.create(user=user,
                             app_name=app_name,
                             name=name,
                             date=timezone.now(),
                             details=details,
                             filename=filename,
                             severity=severity)

    if log:
        severity_display = dict(LogLevel.CHOICES).get(severity, "INFO")
        level = logging.getLevelName(severity_display)
        logging.log(level, e)

    return e


user_logged_in.connect(create_login_event)

