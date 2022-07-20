import itertools
from abc import ABC, abstractmethod
from dataclasses import dataclass
from datetime import timedelta, datetime
from typing import List, Union, Optional

import django.dispatch

from library.log_utils import NotificationBuilder, report_message
from library.utils import segment, flatten_nested_lists


class HealthCheckStat(ABC):
    pass

    def __lt__(self, other):
        return False

    @classmethod
    def to_lines(cls, items: List) -> List[str]:
        return [str(item) for item in items]

    @classmethod
    def sort_order(cls):
        return 0


@dataclass
class HealthCheckActivity(HealthCheckStat):
    emoji: str
    name: str
    amount: int

    def __str__(self):
        amount_str = self.amount
        if self.amount:
            amount_str = f"*{amount_str}*"
        return f"{self.emoji} {self.amount_str} : {self.name}"

    @classmethod
    def sort_order(cls):
        return 1


@dataclass
class HealthCheckAge(HealthCheckStat):
    name: str
    last_performed: Optional[datetime]
    warning_age: timedelta

    # @classmethod
    # def to_lines(cls, items: List) -> List[str]:
    #
    #     items = sorted(items, key=lambda x: (x.warning_age, x.last_performed))
    #     for warning_age, warning_age_group in itertools.groupby(items, lambda x: x.warning_age):

    def __str__(self):
        if not self.last_performed:
            return f":neutral_face: Never Run {self.name}"
        return f":neutral_face: {(datetime.now() - self.last_performed).days} days old : {self.name}"

    @classmethod
    def sort_order(cls):
        return 2


@dataclass
class HealthCheckPercent(HealthCheckStat):
    name: str
    used: str
    available: str

    def __str__(self):
        return f":disk: {self.used} used, {self.available} available : {self.name}"

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


def populate_health_check(notification: NotificationBuilder):
    # calling send_robust will respond with a list of tuples of receiver, the result of calling the receiver
    # the results could be an Exception, a HealthCheckStat, or a list of HealthCheckStats
    # flatten the results, report all the exceptions, then present all the HealthChecks
    stats = flatten_nested_lists([hc[1] for hc in health_check_signal.send_robust(sender=None)])
    checks, exceptions = segment(stats, lambda s: isinstance(s, HealthCheckStat))
    for exec in exceptions:
        notification.add_markdown(f"Exception generating health check: {exec}")

    grouped_checks = [(key, list(values)) for key, values in itertools.groupby(checks, lambda check: type(check))]
    grouped_checks = sorted(grouped_checks, key=lambda gc: gc[0].sort_order())
    lines = list()
    for check_type, checks_typed in grouped_checks:
        lines.extend(check_type.to_lines(checks_typed))

    if lines:
        notification.add_markdown("\n".join(lines), indented=True)