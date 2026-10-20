"""
The event log: Event (severity, app, name, details - what `vg status` lists as recent errors),
ViewEvent (one row per page view), IntegrationActivity (a rolling last-seen per external system,
updated in place through `IntegrationActivity.track`) and ApiRequestLog (one row per call to an
INTEGRATION_API_TRACKING endpoint, written by IntegrationApiMiddleware). Write events with create_event /
create_login_event; the logging handler and middleware in this app fill the rest. An API view links its
ApiRequestLog row to the objects it worked on with log_api_request_objects; what the call wrote is
collected by library/django_utils/change_tracking.py, and any code the call runs can add a line saying
what it did with api_request_note.
"""
import inspect
import json
import logging
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Optional, Union

from django.apps import apps
from django.contrib.auth.models import Group, User
from django.contrib.auth.signals import user_logged_in
from django.contrib.contenttypes.models import ContentType
from django.contrib.postgres.fields import ArrayField
from django.core.serializers.json import DjangoJSONEncoder
from django.db import models
from django.db.models import F, Q
from django.db.models.deletion import SET_NULL
from django.utils import timezone
from model_utils.models import TimeStampedModel

from library.django_utils.change_tracking import active_collector
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


API_KEY_PREFIX = "api:"  # IntegrationActivity.key of an API endpoint, before its view name


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
    # "seen it" on Server Status - errors after this timestamp show again
    last_error_acknowledged = models.DateTimeField(null=True, blank=True)
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

    @property
    def error_acknowledged(self) -> bool:
        if not self.last_error:
            return False
        return bool(self.last_error_acknowledged) and self.last_error_acknowledged >= self.last_error

    @classmethod
    def acknowledge_error(cls, key: str):
        cls.objects.filter(key=key).update(last_error_acknowledged=timezone.now())

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


class ApiRequestLog(models.Model):
    """ One row per call to an INTEGRATION_API_TRACKING endpoint - what was sent and what came back.
        Kept forever unless settings.API_REQUEST_LOG_DAYS is set, then pruned by eventlog.tasks.prune_api_request_log """
    created = models.DateTimeField(auto_now_add=True, db_index=True)
    user = models.ForeignKey(User, null=True, on_delete=SET_NULL)
    method = models.CharField(max_length=10)
    path = models.TextField()
    query_string = models.TextField(blank=True)
    view_name = models.TextField()
    view_kwargs = models.JSONField(default=dict, encoder=DjangoJSONEncoder)  # URL converters can give a UUID
    status_code = models.IntegerField()
    duration_ms = models.IntegerField()
    request_body = models.TextField(null=True)  # trimmed, secrets redacted
    response_body = models.TextField(null=True)  # trimmed, only kept for errors
    # the view's api_log_model, and the pks of those it worked on - empty when it failed before finding one
    content_type = models.ForeignKey(ContentType, null=True, on_delete=SET_NULL)
    object_pks = ArrayField(models.TextField(), default=list)
    # {"app_label.model": {"created": n, "updated": n, "unchanged": n, "deleted": n, "fields": [names]}}
    # - fields are the ones found different on an update, sorted, across the request. {} on rows from before #1413
    changes = models.JSONField(default=dict)
    # created + updated + deleted across all models - the grid's sortable "Changes" column and its "wrote data" filter
    num_changes = models.IntegerField(default=0)
    # what the call did beyond the rows it wrote, from api_request_note - eg data moved onto a re-sent sample sheet
    notes = ArrayField(models.TextField(), default=list)

    class Meta:
        indexes = [models.Index(fields=["view_name", "created"])]

    def __str__(self):
        return f"{self.method} {self.path} ({self.status_code})"

    @property
    def is_error(self) -> bool:
        return self.status_code >= 400

    @property
    def change_rows(self) -> list[dict]:
        """ changes as one dict per model, with its verbose name, for the detail page """
        rows = []
        for label, model_changes in self.changes.items():
            try:
                model_name = apps.get_model(label)._meta.verbose_name
            except (LookupError, ValueError):
                model_name = label  # model since removed
            rows.append({"label": label, "model_name": model_name, **model_changes})
        return rows

    @staticmethod
    def changes_summary(changes: dict) -> str:
        """ One line per model with the counts it has, eg "seqauto.samplesheet: 1 updated (path), 2 unchanged" """
        lines = []
        for label, model_changes in changes.items():
            counts = []
            for outcome in ("created", "updated", "unchanged", "deleted"):
                if count := model_changes.get(outcome):
                    text = f"{count} {outcome}"
                    if outcome == "updated" and (fields := model_changes.get("fields")):
                        text += f" ({', '.join(fields)})"
                    counts.append(text)
            lines.append(f"{label}: {', '.join(counts)}")
        return "\n".join(lines)

    @staticmethod
    def for_object(content_type: ContentType, object_pk: str):
        """ The calls about one object: those that found it, and those for its model that failed before
            finding anything but name its pk in the request body """
        # JSON bodies are stored re-serialised (redacted), so a string appears exactly as json.dumps writes it
        failed_mentioning = Q(object_pks=[], request_body__contains=json.dumps(object_pk))
        return ApiRequestLog.objects.filter(Q(object_pks__contains=[object_pk]) | failed_mentioning,
                                            content_type=content_type)

    @staticmethod
    def has_any_for_model(model: type[models.Model]) -> bool:
        return ApiRequestLog.objects.filter(content_type=ContentType.objects.get_for_model(model)).exists()


def log_api_request_objects(request, *object_pks):
    """ Links this request's ApiRequestLog row to these pks of its view's api_log_model. A DRF Request
        wraps the Django one, which is what IntegrationApiMiddleware reads them from """
    django_request = getattr(request, "_request", request)
    logged_pks = django_request.__dict__.setdefault("api_log_object_pks", set())
    logged_pks.update(str(pk) for pk in object_pks if pk is not None)


def api_request_note(text: str):
    """ Adds a line to the ApiRequestLog row of the tracked API call this runs inside. Code with no request in
        hand (a signal handler) can call it too, as the call's notes live with its change collector """
    logging.info("API request note: %s", text)
    if collector := active_collector():
        collector.notes.append(text)


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
        level = logging.getLevelNamesMapping()[severity_display]
        logging.log(level, e)

    return e


user_logged_in.connect(create_login_event)

