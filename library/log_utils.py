import json
from logging import StreamHandler
import logging
from typing import Dict, Optional

import requests
import rollbar
import sys
import traceback

from django.conf import settings
from django.utils import timezone

from django.contrib.auth.models import User
from rest_framework.request import Request

from eventlog.models import Event
from library.enums.log_level import LogLevel
from threadlocals.threadlocals import get_current_request


def report_event(name: str, request: Request = None, extra_data: Dict = None):
    rollbar.report_message(message=name,
                           level='info',
                           request=request,
                           extra_data=extra_data)

    user: User = None
    details = None
    if request:
        user = request.user

    if extra_data:
        details = extra_data
        if isinstance(extra_data, dict):
            try:
                details = json.dumps(extra_data)
            except:
                details = "(error saving extra_data)"
    elif request:
        details = json.dumps(request.query_params.dict())

    Event.objects.create(user=user,
                         app_name='event',  # could potentially look at the stack trace
                         name=name,
                         date=timezone.now(),
                         details=details,
                         severity=LogLevel.INFO)


def report_message(message: str, level: str = 'warning', request=None, extra_data: dict = None):
    """
    For reporting non-fatal messages to whatever error logging system we want to use. Currently rollbar.
    Note that fatal errors should already be logged with exception catchers elsewhere.
    @param message The message to report on
    @param level The error level, error, warning, info
    @param request the web request (if available)
    @param extra_data a JSON-isable dictionary of extra information
    @param persist_name Should this message be kept permanently, if so give it a name
    """
    print(message)
    if not request:
        request = get_current_request()
    rollbar.report_message(message=message,
                           level=level,
                           request=request,
                           extra_data=extra_data)


def report_exc_info(extra_data=None, request=None):
    if not request:
        request = get_current_request()
    rollbar.report_exc_info(extra_data=extra_data, request=request)
    exc_info = sys.exc_info()
    if exc_info:
        print(exc_info)


def send_notification(message: str, blocks: Optional[Dict] = None, username: Optional[str] = None, emoji: str = ":dna:"):
    """
    Sends a message to your notification service, currently Slack centric.
    If Slack is not configured, this will do nothing.
    @param message The message to send, (if also sending blocks just have message as a summary, wont be displayed)
    @param blocks See https://api.slack.com/messaging/webhooks#advanced_message_formatting
    @param username The username that will appear in Slack
    @param emoji The emoji that will
    """
    sent = False
    if slack := settings.SLACK:
        if slack.get('enabled'):
            if admin_callback_url := slack.get('admin_callback_url'):
                data = {
                    "username": settings.SITE_NAME + (f" {username}" if username else ""),
                    "text": message,
                    "icon_emoji": emoji
                }
                if blocks:
                    data["blocks"] = blocks
                data_str = json.dumps(data)
                print(data_str)
                r = requests.request(
                    headers={"Content-Type": "application/json"},
                    method="POST",
                    url=admin_callback_url,
                    json=data
                )
                try:
                    r.raise_for_status()
                    sent = True
                except:
                    report_exc_info()
    if not sent:
        print("Slack not enabled, did not send message")
        print(message)


def console_logger():
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    handler = StreamHandler()
    handler.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
    logger.addHandler(handler)
    return logger


def get_traceback():
    exec_type, exec_value, _ = sys.exc_info()
    return "".join(traceback.format_exception(exec_value.__class__, exec_value, exec_value.__traceback__))


def log_traceback(level=logging.ERROR):
    tb = get_traceback()
    logging.log(level, tb)
