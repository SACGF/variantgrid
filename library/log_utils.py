import json
import logging
import re
import socket
import sys
import traceback
from abc import ABC, abstractmethod
from logging import StreamHandler
from re import Match
from typing import Optional, Any, Union

import markdown
import requests
import rollbar
from django.conf import settings
from django.contrib.admin.options import get_content_type_for_model
from django.contrib.auth.models import User
from django.db.models import Model
from django.forms import ModelForm
from django.utils import timezone
from markdown import markdown
from requests import HTTPError
from rest_framework.request import Request
from threadlocals.threadlocals import get_current_request, get_current_user

from eventlog.models import Event
from library.constants import MINUTE_SECS
from library.enums.log_level import LogLevel
from library.utils import pretty_label


def get_current_logged_in_user() -> Optional[User]:
    user: User = get_current_user()
    if user:
        if not user.is_authenticated:
            return None
    return user


def report_event(name: str, request: Request = None, extra_data: dict = None):
    rollbar.report_message(message=name,
                           level='info',
                           request=request,
                           extra_data=extra_data)

    if request is None:
        request = get_current_request()

    user: Optional[User] = None
    details = None
    if request and request.user.is_authenticated:
        user = request.user

    if extra_data:
        details = extra_data
        if isinstance(extra_data, dict):
            try:
                details = json.dumps(extra_data)
            except TypeError:
                details = "(error saving extra_data)"
    elif request:
        details = json.dumps({**request.POST.dict(), **request.GET.dict()})

    Event.objects.create(user=user,
                         app_name='event',  # could potentially look at the stack trace
                         name=name,
                         date=timezone.now(),
                         details=details,
                         severity=LogLevel.INFO)


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
    from django.contrib.admin.models import CHANGE, LogEntry
    if not user:
        user = get_current_logged_in_user()

    if isinstance(message, dict):
        message = json.dumps(message)

    return LogEntry.objects.log_action(
        user_id=user.pk,
        content_type_id=get_content_type_for_model(obj).pk,
        object_id=obj.pk,
        object_repr=str(obj),
        action_flag=CHANGE,
        change_message=message,
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
    print(exception_message)
    extra_data["exception_message"] = exception_message

    if not request:
        request = get_current_request()

    rollbar.report_message(message=message,
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


class NotificationBuilder:
    # Preference is to use AdminNotificationBuilder or LabNotificationBuilder now

    SLACK_EMOJI_RE = re.compile(r":[A-Z_]+:", re.IGNORECASE)
    DE_P = re.compile(r"<p>(.*?)</p>", re.IGNORECASE | re.DOTALL)
    LINK_RE = re.compile(r"<(?P<link>http.*)\|(?P<name>.*)>")

    @staticmethod
    def slack_markdown_to_html(markdown_txt: str, surround_with_div: bool = False):
        if markdown_txt:
            if not isinstance(markdown_txt, str):
                raise TypeError(f"Expected string or bytes-like object, got {markdown_txt}")
            markdown_txt = re.sub(NotificationBuilder.SLACK_EMOJI_RE, "", markdown_txt)

            # convert Slack links <http://foobar.com|FooBar> to proper Markdown links [FooBar](http://foobar.com)
            def link_fixer(match: Match):
                return f"[{match.group('name')}]({match.group('link')})"

            markdown_txt = re.sub(NotificationBuilder.LINK_RE, link_fixer, markdown_txt)

            markdown_html = markdown(markdown_txt)
            if surround_with_div:
                markdown_html = f"<div>{markdown_html}</div>"
            else:
                if de_ped := NotificationBuilder.DE_P.match(markdown_html):
                    markdown_html = de_ped.group(1)

            markdown_html = markdown_html.replace("\n", "<br/>")
            return markdown_html
        return ""

    # Different kind of blocks

    class Block(ABC):

        @abstractmethod
        def as_text(self) -> str:
            raise NotImplementedError(f"{self} has not implemented as_text method")

        @abstractmethod
        def as_html(self) -> str:
            raise NotImplementedError(f"{self} has not implemented as_html method")

        @abstractmethod
        def as_slack(self) -> Union[dict, list]:
            raise NotImplementedError(f"{self} has not implemented as_slack method")

    class HeaderBlock(Block):

        def __init__(self, header_text: str):
            self.header_text = header_text

        def as_text(self):
            return self.header_text

        def as_html(self):
            return f"<h4>{NotificationBuilder.slack_markdown_to_html(self.header_text)}</h4>"

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
            self.fields: list[Union[str, Any]] = []

        def add_field(self, label: str, value: str):
            self.fields.append((label, value))

        def as_text(self):
            def as_line(field: tuple[str, Any]):
                return f"{field[0]}: {field[1]}"
            return "\n".join(as_line(field) for field in self.fields)

        def as_html(self):
            def as_field(field: tuple[str, Any]):
                return f"<p><span style='color:#888'>{NotificationBuilder.slack_markdown_to_html(field[0])}:</span><br/><span style='font-family:monospace'>{NotificationBuilder.slack_markdown_to_html(field[1])}</span></p>"

            return "".join([as_field(field) for field in self.fields])

        def as_slack(self):
            blocks = []

            def as_field(the_field: tuple[str, Any]):
                return {
                    "type": "mrkdwn",
                    "text": f"*{the_field[0]}:*\n{the_field[1]}"
                }

            def add_field_list(fields: list[tuple[str, Any]]):
                if fields:
                    blocks.append({
                        "type": "section",
                        "fields": [as_field(the_field) for the_field in fields]
                    })

            field_list = []
            for field in self.fields:
                field_list.append(field)
                if len(field_list) == 2:
                    add_field_list(field_list)
                    field_list = []
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

    def __init__(self, message: str):
        """
        :param message: WARNING, this is ignored for Slack Notifications (maybe turn into a header?)
        """
        self.message = message
        self.blocks: list[NotificationBuilder.Block] = []
        self.sent = False

    def _last_block(self) -> Optional[Block]:
        if self.blocks:
            return self.blocks[-1]
        return None

    def add_header(self, text) -> 'NotificationBuilder':
        self.blocks.append(NotificationBuilder.HeaderBlock(text))
        return self

    def add_field(self, label: str, value: Union[str, int]) -> 'NotificationBuilder':
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

    def merge(self, *bulk_notifications: 'NotificationBuilder') -> 'NotificationBuilder':
        for bulk_notification in bulk_notifications:
            if bulk_notification == self:
                raise ValueError("Can't merge NotificationBuilder with itself")
            bulk_notification.sent = True
            self.blocks.extend(bulk_notification.blocks)

            # Not sure what the logic of the below was, only include a header if it's the first row?

            # if bulk_notification.blocks and isinstance(bulk_notification.blocks[0], NotificationBuilder.HeaderBlock):
            #     self.blocks.append(bulk_notification.blocks[0])
            #
            # for block in bulk_notification.blocks:
            #     if isinstance(block, (NotificationBuilder.FieldsBlock, NotificationBuilder.MarkdownBlock)):
            #         self.blocks.append(block)
        return self

    @property
    def webhook_url(self) -> Optional[str]:
        return None

    def as_slack(self) -> list:
        slack_blocks: list = []
        for block in self.blocks:
            slack_bit = block.as_slack()
            if isinstance(slack_bit, dict):
                slack_blocks.append(slack_bit)
            elif isinstance(slack_bit, list):
                slack_blocks += slack_bit
        return slack_blocks

    def as_html(self):
        return "<br/>".join([block.as_html() for block in self.blocks])

    def as_text(self):
        return "\n".join([block.as_text() for block in self.blocks])

    def send(self):
        self.sent = True
        Event.objects.create(user=get_current_logged_in_user(),
                             app_name='event',  # could potentially look at the stack trace
                             name=self.message,
                             date=timezone.now(),
                             details=self.as_text(),
                             severity=LogLevel.INFO)
        send_notification(message=self.message, blocks=self.as_slack(), slack_webhook_url=self.webhook_url)

    def __del__(self):
        if not self.sent:
            report_message(f"Created a NotificationBuilder but did not call send {self.message}")


class AdminNotificationBuilder(NotificationBuilder):

    def __init__(self, message: str, is_communication: bool = False):
        self.is_communication = is_communication
        super().__init__(message=message)

    def send(self):
        from email_manager.models import EmailLog
        super().send()
        if self.is_communication and settings.ADMIN_EMAIL_NOTIFICATION:
            email_recipients = [user.email for user in User.objects.filter(is_superuser=True) if user.email]
            EmailLog.send_mail(subject=self.message,
                               html=self.as_html(),
                               text=self.as_text(),
                               from_email=settings.ADMIN_EMAIL_NOTIFICATION,
                               recipient_list=email_recipients,
                               allow_users_to_see_others=True)


def slack_bot_username():
    env_name = socket.gethostname().lower().split('.')[0].replace('-', '')
    site_name = settings.SITE_NAME
    if site_name.replace(' ', '').lower() == env_name:
        return site_name
    else:
        return f"{site_name} ({env_name})"


# message limit is actually 4000, but this gives us a lot of leway
SLACK_CHARACTER_LIMIT = 3500


def send_notification(
        message: str,
        blocks: Optional[list[NotificationBuilder.Block]] = None,
        slack_webhook_url: Optional[str] = None):
    """
    Sends a message to your notification service, currently Slack centric.
    Best practise is to use a NotificationBuilder with send() rather than calling send_notification directly

    If Slack is not configured, this will do nothing.
    :param message: appears to be IGNORED, should probably do something about that
    :param blocks: See https://api.slack.com/messaging/webhooks#advanced_message_formatting
    :param slack_webhook_url: Provide the slack URL, if not provided will get from settings (if enabled)
    """
    emoji = ':dna:'
    if slack := settings.SLACK:
        if slack_webhook_url is None:
            if slack.get('enabled'):
                slack_webhook_url = slack.get('admin_callback_url')
        if env_emoji := slack.get('emoji'):
            emoji = env_emoji

    if slack_webhook_url:
        username = slack_bot_username()

        data = {
            "username": username,
            "text": message,
            "icon_emoji": emoji
        }
        character_count = 0
        data_blocks = []
        for block in blocks:
            this_block = json.dumps(block)
            character_count += len(this_block)
            if character_count >= SLACK_CHARACTER_LIMIT:
                data_blocks.append(
                    NotificationBuilder.MarkdownBlock("Message is too big for slack, check EventLog for the full message").as_slack()
                )
                break
            data_blocks.append(block)

        if data_blocks:
            data["blocks"] = data_blocks

        r = requests.post(
            headers={"Content-Type": "application/json"},
            url=slack_webhook_url,
            json=data,
            timeout=MINUTE_SECS,
        )
        try:
            r.raise_for_status()
        except HTTPError:
            # give Rollbar the Slack JSON if Slack doesn't want it
            report_exc_info(extra_data=data)
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
    LogLevel.DEBUG: 'info',
    LogLevel.INFO: 'info',
    LogLevel.WARNING: 'warning',
    LogLevel.ERROR: 'danger'
}


def log_level_to_int(log_level: LogLevel) -> int:
    return _LOG_LEVEL_SEVERITY.get(log_level, 0)


def log_level_to_bootstrap(log_level: LogLevel):
    return _LOG_LEVEL_TO_BOOTSTRAP[log_level]
