import json
from logging import StreamHandler
import logging
import socket
import markdown
from typing import Dict, Optional, List, Tuple, Any, Union

import requests
import rollbar
import sys
import traceback
import re

from django.conf import settings
from django.utils import timezone

from django.contrib.auth.models import User
from markdown import markdown
from rest_framework.request import Request

from eventlog.models import Event
from library.enums.log_level import LogLevel
from threadlocals.threadlocals import get_current_request, get_current_user


def report_event(name: str, request: Request = None, extra_data: Dict = None):
    rollbar.report_message(message=name,
                           level='info',
                           request=request,
                           extra_data=extra_data)

    if request is None:
        request = get_current_request()

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

    SLACK_EMOJI_RE = re.compile(r"[:][A-Z_]+[:]", re.IGNORECASE)

    @staticmethod
    def slack_markdown_to_html(markdown_txt: str, surround_with_div: bool = False):
        if markdown_txt:
            markdown_txt = re.sub(NotificationBuilder.SLACK_EMOJI_RE, "", markdown_txt)
            markdown_html = markdown(markdown_txt)
            if surround_with_div:
                markdown_html = f"<div>{markdown_html}</div>"
            return markdown_html
        return ""

    # Different kind of blocks

    class Block:

        def as_text(self) -> str:
            raise NotImplementedError(f"{self} has not implemented as_text method")

        def as_html(self) -> str:
            raise NotImplementedError(f"{self} has not implemented as_html method")

        def as_slack(self) -> Union[Dict, List]:
            raise NotImplementedError(f"{self} has not implemented as_slack method")

    class HeaderBlock(Block):

        def __init__(self, header_text: str):
            self.header_text = header_text

        def as_text(self):
            return self.header_text

        def as_html(self):
            return f"<h4>{self.header_text}</h4>"

        def as_slack(self):
            return {
                "type": "header",
                "text": {
                    "type": "plain_text",
                    "text": self.header_text,
                    "emoji": True
                }
            }

    class FieldsBlock(Block):

        def __init__(self):
            self.fields: List[Union[str, Any]] = list()

        def add_field(self, label: str, value: str):
            self.fields.append((label, value))

        def as_text(self):
            def as_line(field: Tuple[str, Any]):
                return f"{field[0]}: {field[1]}"
            return "\n".join(as_line(field) for field in self.fields)

        def as_html(self):
            def as_field(field: Tuple[str, Any]):
                return f"<label>{NotificationBuilder.slack_markdown_to_html(field[0])}</label>: {field[1]}"
            return "<div>" + "<br/>".join([as_field(field) for field in self.fields]) + "</div>"

        def as_slack(self):
            blocks = list()

            def as_field(field: Tuple[str, Any]):
                return {
                    "type": "mrkdwn",
                    "text": f"*{field[0]}:*\n{field[1]}"
                }

            def add_field_list(fields: List[Tuple[str, Any]]):
                if fields:
                    blocks.append({
                        "type": "section",
                        "fields": [as_field(field) for field in fields]
                    })

            field_list = list()
            for field in self.fields:
                field_list.append(field)
                if len(field_list) == 2:
                    add_field_list(field_list)
                    field_list = list()
            add_field_list(field_list)
            return blocks

    class MarkdownBlock(Block):

        def __init__(self, markdown_txt: str, indented: bool = False):
            self.markdown_txt = markdown_txt
            self.indented = indented

        def as_text(self):
            text = self.markdown_txt
            if self.indented:
                text = "    " + (text or "<No Data>")
            return text

        def as_slack(self):
            text = self.markdown_txt
            if self.indented:
                text = ">>> " + (text or "_No Data_")
            return {"type": "section", "text": {"type": "mrkdwn", "text": text}}

        def as_html(self):
            return NotificationBuilder.slack_markdown_to_html(self.markdown_txt)

    class DividerBlock(Block):

        def __init__(self):
            pass

        def as_text(self):
            return "--------------"

        def as_slack(self):
            return {"type": "divider"}

        def as_html(self):
            return "<hr/>"

    # end blocks

    def __init__(self, message: str, emoji: str = ":dna:"):
        self.message = message
        self.emoji = emoji
        self.blocks: List[NotificationBuilder.Block] = list()
        self.sent = False

    def _last_block(self) -> Optional[Block]:
        if self.blocks:
            return self.blocks[-1]
        return None

    def add_header(self, text) -> 'NotificationBuilder':
        self.blocks.append(NotificationBuilder.HeaderBlock(text))
        return self

    def add_field(self, label: str, value: str) -> 'NotificationBuilder':
        fields_block: Optional[NotificationBuilder.FieldsBlock] = None
        if last_block := self._last_block():
            if isinstance(last_block, NotificationBuilder.FieldsBlock):
                fields_block = last_block
        if not fields_block:
            fields_block = NotificationBuilder.FieldsBlock()
            self.blocks.append(fields_block)
        fields_block.add_field(label=label, value=value)
        return self

    def add_divider(self) -> 'NotificationBuilder':
        self.blocks.append(NotificationBuilder.DividerBlock())
        return self

    def add_markdown(self, text, indented=False):
        self.blocks.append(NotificationBuilder.MarkdownBlock(markdown_txt=text, indented=indented))
        return self

    @property
    def webhook_url(self) -> Optional[str]:
        return None

    def as_slack(self):
        slack_blocks: List = list()
        for block in self.blocks:
            slack_bit = block.as_slack()
            if isinstance(slack_bit, dict):
                slack_blocks.append(slack_bit)
            elif isinstance(slack_bit, list):
                slack_blocks += slack_bit

        env_name = socket.gethostname().lower().split('.')[0].replace('-', '')
        slack_blocks.append({
            "type": "context",
            "elements": [
                {
                    "type": "plain_text",
                    "text": env_name,
                    "emoji": False
                }
            ]
        })
        return slack_blocks

    def as_html(self):
        return "<br/>".join([block.as_html() for block in self.blocks])

    def as_text(self):
        return "\n".join([block.as_text() for block in self.blocks])

    def send(self):
        self.sent = True
        Event.objects.create(user=get_current_user(),
                             app_name='event',  # could potentially look at the stack trace
                             name=self.message,
                             date=timezone.now(),
                             details=self.as_text(),
                             severity=LogLevel.INFO)
        send_notification(message=self.message, blocks=self.as_slack(), emoji=self.emoji, slack_webhook_url=self.webhook_url)

    def __del__(self):
        if not self.sent:
            report_message(f"Created a NotificationBuilder but did not call send {self.message}")


def send_notification(
        message: str,
        blocks: Optional[Dict] = None,
        username: Optional[str] = None,
        emoji: str = ":dna:",
        slack_webhook_url: Optional[str] = None):
    """
    Sends a message to your notification service, currently Slack centric.
    Best practise is to use a NotificationBuilder with send() rather than calling send_notification directly

    If Slack is not configured, this will do nothing.
    @param message The message to send, (if also sending blocks just have message as a summary, wont be displayed)
    @param blocks See https://api.slack.com/messaging/webhooks#advanced_message_formatting
    @param username The username that will appear in Slack
    @param emoji The emoji that will
    @param slack_webhook_url Provide the slack URL, if not provided will get from settings (if enabled)
    """
    sent = False
    if slack_webhook_url is None:
        if slack := settings.SLACK:
            if slack.get('enabled'):
                slack_webhook_url = slack.get('admin_callback_url')

    if slack_webhook_url:
        data = {
            "username": (settings.SITE_NAME + (f" {username}" if username else "")).strip(),
            "text": message,
            "icon_emoji": emoji
        }
        if blocks:
            data["blocks"] = blocks
        data_str = json.dumps(data)

        r = requests.request(
            headers={"Content-Type": "application/json"},
            method="POST",
            url=slack_webhook_url,
            json=data
        )
        try:
            r.raise_for_status()
            sent = True
        except:
            report_exc_info()
    else:
        # fallback to Rollbar if Slack isn't configured
        report_event(name=message)


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
