"""
Telling Rollbar and the logs what happened: report_exc_info / report_message for caught problems,
log_admin_change / log_saved_form for audit, and the console and traceback helpers. Messages below
settings.ROLLBAR['min_level'] are not sent to Rollbar. The current request is found through
library/request_context.py. Notifying humans (Slack, email, EventLog) is eventlog/notifications.py.
"""
import json
import logging
import sys
import traceback
from logging import StreamHandler
from typing import Optional, Union

import rollbar
from django.conf import settings
from django.contrib.admin.models import CHANGE, LogEntry
from django.contrib.auth.models import User
from django.db.models import Model
from django.forms import ModelForm

from library.enums.log_level import LogLevel
from library.request_context import get_current_request, get_current_user
from library.utils import pretty_label


def get_current_logged_in_user() -> Optional[User]:
    user: User = get_current_user()
    if user:
        if not user.is_authenticated:
            return None
    return user


_ROLLBAR_LEVEL_SEVERITY = {level: i for i, level in enumerate(['debug', 'info', 'warning', 'error', 'critical'])}


def rollbar_report_message(message: str, level: str, request=None, extra_data: Optional[dict] = None):
    """ Sends to Rollbar unless level is below settings.ROLLBAR['min_level'] """
    if min_level := settings.ROLLBAR.get('min_level'):
        max_severity = len(_ROLLBAR_LEVEL_SEVERITY)
        if _ROLLBAR_LEVEL_SEVERITY.get(level, max_severity) < _ROLLBAR_LEVEL_SEVERITY.get(min_level, 0):
            return
    rollbar.report_message(message=message,
                           level=level,
                           request=request,
                           extra_data=extra_data)




def log_saved_form(form: ModelForm, user: Optional[User] = None):
    if form.has_changed() and form.is_valid():
        changed_fields_str = ", ".join([pretty_label(f) for f in sorted(form.changed_data)])
        message = f"Changed {changed_fields_str} via form."
        log_admin_change(obj=form.instance, message=message, user=user)


def log_admin_change(obj: Model, message: Union[str, dict], user: Optional[User] = None):
    """
    Log that an object has been successfully changed.

    The default implementation creates an admin LogEntry object.
    """
    if not user:
        user = get_current_logged_in_user()

    if isinstance(message, dict):
        message = json.dumps(message)

    return LogEntry.objects.log_actions(
        user_id=user.pk,
        queryset=[obj],
        #content_type_id=get_content_type_for_model(obj).pk,
        #object_id=obj.pk,
        #object_repr=str(obj),
        action_flag=CHANGE,
        change_message=message,
        single_object=True
    )


def report_message(message: str, level: str = 'warning', request=None, extra_data: dict = None):
    """
    For reporting non-fatal messages to whatever error logging system we want to use. Currently rollbar.
    Note that fatal errors should already be logged with exception catchers elsewhere.
    @param message The message to report on
    @param level The error level, error, warning, info
    @param request the web request (if available)
    @param extra_data a JSON-able dictionary of extra information
    """
    if not extra_data:
        extra_data = {}
    target = extra_data.get("target")
    exception_message = message
    if target:
        exception_message += f": {target}"
    logging.info(exception_message)
    extra_data["exception_message"] = exception_message

    if not request:
        request = get_current_request()

    rollbar_report_message(message=message,
                            level=level,
                            request=request,
                            extra_data=extra_data)


def report_exc_info(extra_data=None, request=None, report_externally=True):
    if not request:
        request = get_current_request()

    exc_info = sys.exc_info()
    if report_externally:
        if not extra_data:
            extra_data = {}
        extra_data["exception_message"] = str(exc_info[1])
        rollbar.report_exc_info(exc_info=exc_info, extra_data=extra_data, request=request)

    if exc_info:
        traceback.print_exc()




def console_logger():
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    handler = StreamHandler()
    handler.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
    logger.addHandler(handler)
    return logger


def get_traceback():
    _, exec_value, _ = sys.exc_info()
    return "".join(traceback.format_exception(exec_value.__class__, exec_value, exec_value.__traceback__))


def log_traceback(level=logging.ERROR):
    tb = get_traceback()
    logging.log(level, tb)


_LOG_LEVEL_SEVERITY = {
    LogLevel.DEBUG: 0,
    LogLevel.INFO: 1,
    LogLevel.WARNING: 2,
    LogLevel.ERROR: 3
}

_LOG_LEVEL_TO_BOOTSTRAP = {
    None: 'ok',
    LogLevel.DEBUG: 'info',
    LogLevel.INFO: 'info',
    LogLevel.WARNING: 'warning',
    LogLevel.ERROR: 'danger'
}


def log_level_to_int(log_level: LogLevel) -> int:
    return _LOG_LEVEL_SEVERITY.get(log_level, 0)


def log_level_to_bootstrap(log_level: LogLevel):
    return _LOG_LEVEL_TO_BOOTSTRAP[log_level]
