import itertools
from abc import ABC, abstractmethod
from dataclasses import dataclass
from datetime import timedelta, datetime
from typing import List, Union, Optional

import django.dispatch
from django.db.models import Model, Q
from django.utils.timezone import localtime

from library.log_utils import NotificationBuilder, report_message
from library.utils import segment, flatten_nested_lists, model_has_field


class HealthCheckRequest:

    def __init__(self, since: datetime, now: Optional[datetime] = None):
        self.since = since
        self.now = now or datetime.now()


class HealthCheckStat(ABC):

    def __lt__(self, other):
        return False

    @classmethod
    def to_lines(cls, items: List, health_request: HealthCheckRequest) -> List[str]:
        return [str(item) for item in items]

    @classmethod
    def sort_order(cls):
        return 0


@dataclass
class HealthCheckRecentActivity(HealthCheckStat):
    emoji: str
    name: str
    amount: int
    sub_type: Optional[str] = None
    extra: Optional[str] = None
    stand_alone: bool = False

    # Consider this taking a QuerySet of the objects instead
    # that way Server Status web page could use them

    def __str__(self):
        amount_str = self.amount
        if self.amount:
            amount_str = f"*{amount_str}*"
        result = f"{self.emoji} {amount_str} : {self.name}"
        if self.sub_type:
            result = f"{result} {self.sub_type}"
        if self.extra:
            result = f"{result} : {self.extra}"
        return result

    @classmethod
    def sort_order(cls):
        return 1

    @staticmethod
    def simple_report(
            health_request: HealthCheckRequest,
            model: Model,
            emoji: str,
            created: bool = False,
            modified: bool = False) -> List['HealthCheckRecentActivity']:
        name = model.__name__
        name += "es" if name.endswith("s") else "s"

        if verbose_name := model._meta.verbose_name_plural:
            if verbose_name[0].isupper():
                # assume verbose name plural is auto generated if it doesn't start with a capital
                name = verbose_name

        response = list()
        if created:
            q: Q
            if model_has_field(model, 'created'):
                q = Q(created__gte=health_request.since) & Q(created__lte=health_request.now)
            else:
                # some older models have "date" rather than being TimestampModel
                # maybe we can migrate them to timestamp model and populate the created field?
                q = Q(date__gte=health_request.since) & Q(date__lte=health_request.now)
            response.append(HealthCheckRecentActivity(
                emoji=emoji,
                name=name,
                amount=model.objects.filter(q).count(),
                sub_type="Created"
            ))
        if modified:
            # don't bother supporting 'date' instead of 'created' here as we need
            # modified as well.
            response.append(HealthCheckRecentActivity(
                emoji=emoji,
                name=name,
                amount=model.objects.filter(
                    created__lt=health_request.since,
                    modified__gte=health_request.since,
                    modified__lte=health_request.now).count(),

                sub_type="Modified"
            ))
        return response


@dataclass
class HealthCheckTotalAmount(HealthCheckStat):
    emoji: str
    name: str
    amount: int
    extra: Optional[str] = None

    def __str__(self):
        amount_str = self.amount
        if self.amount:
            amount_str = f"*{amount_str}*"
        result = f"{self.emoji} {amount_str} : {self.name}"
        if self.extra:
            result = f"{result} - {self.extra}"
        return result

    @classmethod
    def sort_order(cls):
        return 1


@dataclass
class HealthCheckAge(HealthCheckStat):
    name: str
    last_performed: Optional[datetime]
    warning_age: timedelta

    @property
    def last_performed_tz(self):
        if self.last_performed:
            return localtime(self.last_performed)

    @classmethod
    def to_lines(cls, items: List, health_request: HealthCheckRequest) -> List[str]:
        return [item.str_with_now(health_request.now) for item in items]

    def str_with_now(self, now: datetime):
        if not self.last_performed_tz:
            return f":dizzy_face: Never Run : {self.name}"
        return f":neutral_face: {(now - self.last_performed_tz).days} days old : {self.name}"

    def __str__(self):
        return self.str_with_now(datetime.now())

    @classmethod
    def sort_order(cls):
        return 2


@dataclass
class HealthCheckCapacity(HealthCheckStat):
    name: str
    used: str
    available: str
    warning: bool = False

    def __str__(self):
        emoji = ":floppy_disk:" if not self.warning else ":fire:"
        return f"{emoji} {self.name} ({self.used} used, {self.available} available)"

    @classmethod
    def sort_order(cls):
        return 3


@dataclass
class HealthCheckCustom(HealthCheckStat):
    text: str

    def __str__(self):
        return self.text

    @classmethod
    def sort_order(cls):
        return 4


health_check_signal = django.dispatch.Signal()


def populate_health_check(notification: NotificationBuilder, since: Optional[datetime] = None):
    # calling send_robust will respond with a list of tuples of receiver, the result of calling the receiver
    # the results could be an Exception, a HealthCheckStat, or a list of HealthCheckStats
    # flatten the results, report all the exceptions, then present all the HealthChecks

    now = localtime()
    if since is None:
        since = now - timedelta(days=1)
    health_request = HealthCheckRequest(since=since, now=now)

    stats = flatten_nested_lists([hc[1] for hc in health_check_signal.send_robust(sender=None, health_request=health_request)])
    checks, exceptions = segment(stats, lambda s: isinstance(s, HealthCheckStat))
    for exec in exceptions:
        notification.add_markdown(f"Exception generating health check: {exec}")

    grouped_checks = [(key, list(values)) for key, values in itertools.groupby(checks, lambda check: type(check))]
    grouped_checks = sorted(grouped_checks, key=lambda gc: gc[0].sort_order())
    lines = list()
    for check_type, checks_typed in grouped_checks:
        lines.extend(check_type.to_lines(checks_typed, health_request=health_request))

    if lines:
        notification.add_markdown("\n".join(lines), indented=True)