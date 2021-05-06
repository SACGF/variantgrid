import json
from logging import StreamHandler
import logging
import socket
from typing import Dict, Optional, List, Tuple, Any

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


class NotificationBuilder:

    def __init__(self, message: str, emoji: str = ":dna:"):
        self.message = message
        self.emoji = emoji
        self.blocks = list()
        self.sent = False
        self.fields: List[Tuple[str, Any]] = list()

    def add_header(self, text) -> 'NotificationBuilder':
        self._check()
        self.blocks.append({
            "type": "header",
            "text": {
                "type": "plain_text",
                "text": text,
                "emoji": True
            }
        })
        return self

    def _check(self):
        self._field_check()

    def _field_check(self):
        def as_field(field: Tuple[str, Any]):
            return {
                "type": "mrkdwn",
                "text": f"*{field[0]}:*\n{field[1]}"
            }

        def add_field_list(fields: List[Tuple[str, Any]]):
            if fields:
                self.blocks.append({
                    "type": "section",
                    "fields": [as_field(field) for field in fields]
                })

        if fields := self.fields:
            field_list = list()
            for field in fields:
                field_list.append(field)
                if len(field_list) == 2:
                    add_field_list(field_list)
                    field_list = list()
            add_field_list(field_list)
            self.fields = list()

    def add_field(self, label: str, value: str, single=False):
        self.fields.append((label, value))

    def add_divider(self) -> 'NotificationBuilder':
        self._check()
        self.blocks.append({"type": "divider"})
        return self

    def add_markdown(self, text, indented=False):
        self._check()
        if indented:
            text = ">>> " + (text or "_No Data_")
        self.blocks.append({"type": "section", "text": {"type": "mrkdwn", "text": text}})
        return self

    def send(self):
        self._check()
        self.sent = True
        if blocks := self.blocks:
            env_name = socket.gethostname().lower().split('.')[0].replace('-', '')
            blocks.append({
                "type": "context",
                "elements": [
                    {
                        "type": "plain_text",
                        "text": env_name,
                        "emoji": False
                    }
                ]
            })
            send_notification(message=self.message, blocks=blocks, emoji=self.emoji)

    def __del__(self):
        if not self.sent:
            report_message(f"Created a NotificationBuilder but did not call send {self.message}")


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
                    "username": (settings.SITE_NAME + (f" {username}" if username else "")).strip(),
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
