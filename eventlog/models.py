import inspect
import logging
from typing import Optional

from django.contrib.auth.models import User
from django.contrib.auth.signals import user_logged_in
from django.db import models
from django.db.models.deletion import SET_NULL
from django.utils import timezone
from model_utils.models import TimeStampedModel

from library.enums.log_level import LogLevel
from library.utils import empty_dict


class ViewEvent(TimeStampedModel):
    user = models.ForeignKey(User, null=True, on_delete=SET_NULL)
    view_name = models.TextField()
    args = models.JSONField(null=False, blank=True, default=empty_dict)
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

    def can_write(self, user):
        return user.is_superuser

    def __str__(self):
        MAX_DETAILS_LENGTH = 200
        details = self.details
        if details and (len(details) > MAX_DETAILS_LENGTH):
            details = details[:MAX_DETAILS_LENGTH // 2] + " ... " + details[-MAX_DETAILS_LENGTH // 2:]
        user_msg = ''
        if self.user:
            user_msg = f"User: {self.user}, "

        return f"{user_msg}{self.app_name}/{self.name}: {details}"


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
        log=True):

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


user_logged_in.connect(create_login_event)
