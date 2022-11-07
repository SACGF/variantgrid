import itertools
from abc import ABC
from dataclasses import dataclass
from datetime import timedelta, datetime
from typing import List, Union, Optional, Type

import django.dispatch
from django.db.models import Model, Q
from django.utils.timezone import localtime

from library.log_utils import NotificationBuilder
from library.utils import flatten_nested_lists, model_has_field

"""
HealthChecks are generated nightly and posted in Slack.
In future these classes might be extended to show the data live on a webpage
"""


class HealthCheckRequest:
    """
    Class that gives context of the HealthRequest, specifically the time window that we care about
    """

    def __init__(self, since: datetime, now: Optional[datetime] = None):
        self.since = since
        self.now = now or datetime.now()


class HealthCheckStat(ABC):

    def __lt__(self, other):
        return False

    @classmethod
    def to_lines(cls, items: List, health_request: HealthCheckRequest) -> List[str]:
        # given a list of items of this class, generate lines of output
        return [str(item) for item in items]

    @classmethod
    def is_recent_activity(cls):
        # indicates if this state represents time within a window (e.g. 3 records updates in 24 hours)
        # or more of a snapshot state (e.g. 10,343 records in total)
        return False

    @classmethod
    def sort_order(cls):
        # how this stat should be sorted with other stats
        return 0


@dataclass
class HealthCheckRecentActivity(HealthCheckStat):
    emoji: str  # The Slack Emoji to put next to the message
    name: str  # The name of the model we're reporting on
    amount: Union[int, str]   # how many records were created/updated/deleted etc
    sub_type: Optional[str] = None   # if we're talking about created/updated or deleted
    extra: Optional[str] = None   # extra text (not used for grouping)
    stand_alone: bool = False  # should this always be reported on in its own line

    # Consider this taking a QuerySet of the objects instead
    # that way Server Status web page could use them

    @classmethod
    def is_recent_activity(cls):
        return True

    @property
    def is_zero(self):
        # if is_zero (and not stand_alone) then the record can be summarised
        return self.amount == 0 or self.amount is None or self.amount == ""

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
            model: Type[Model],
            emoji: str,
            created: bool = False,
            modified: bool = False) -> List['HealthCheckRecentActivity']:
        """
        Makes a recent activity record for a given model
        :param health_request: The request that gives us the time window
        :param model: The model to report create and/or updated
        :param emoji: Emoji that represents the model
        :param created: If true, report on new records : model must have a column called created or date
        :param modified: If true, report on modified records (that weren't also created in the time window) : model must have a column called modified
        :return: One or two health check stats depending on if created and/or modified are true
        """
        name = model.__name__
        name += "es" if name.endswith("s") else "s"

        if verbose_name := model._meta.verbose_name_plural:
            if verbose_name[0].isupper():
                # assume verbose name plural is auto generated if it doesn't start with a capital
                name = verbose_name

        response = []
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

    @classmethod
    def to_lines(cls, items: List['HealthCheckRecentActivity'], health_request: HealthCheckRequest) -> List[str]:
        """
        Display standalone lines at the front, and then non-zero items per line, then all the zeros in one line
        """
        items = sorted(items, key=lambda hc: (hc.name, hc.sub_type))
        zeros: List[HealthCheckRecentActivity] = []
        stand_alone: List[HealthCheckRecentActivity] = []
        values: List[HealthCheckRecentActivity] = []
        for item in items:
            if item.stand_alone:
                stand_alone.append(item)
            elif item.is_zero:
                zeros.append(item)
            else:
                values.append(item)

        output = []
        for stand_alone_hc in stand_alone:
            output.append(str(stand_alone_hc))
        for valued in values:
            output.append(str(valued))
        if zeros:
            flattened_empty = []
            for name, hc_sub_types in itertools.groupby(zeros, key=lambda hc: hc.name):
                flattened_empty.append(name + " " + "/".join(hc.sub_type for hc in hc_sub_types))
            all_empty = ", ".join(flattened_empty)
            output.append(f":open_file_folder: 0 : {all_empty}")
        return output


@dataclass
class HealthCheckTotalAmount(HealthCheckStat):
    # Represents a total amount, e.g. number of records in total
    emoji: str
    name: str
    amount: int
    extra: Optional[str] = None

    def __str__(self):
        amount_str = self.amount
        if self.amount:
            amount_str = f"*{amount_str:,}*"
        result = f"{self.emoji} {amount_str} : {self.name}"
        if self.extra:
            result = f"{result} - {self.extra}"
        return result

    @classmethod
    def sort_order(cls):
        return 2


@dataclass
class HealthCheckCapacity(HealthCheckStat):
    # Used to represent available disk space
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
class HealthCheckAge(HealthCheckStat):
    # Used to report on datasets that should be periodically refreshed
    name: str
    now: datetime
    last_performed: Optional[datetime]  # when was this last run, None if never run
    warning_age: timedelta  # how old can this item be before the user should be warned, causes increasingly angry faces when exceeded

    _NEVER_RUN_AGE = 100000

    @property
    def last_performed_tz(self):
        if self.last_performed:
            return localtime(self.last_performed)

    @property
    def age(self) -> Optional[timedelta]:
        if last_performed_tz := self.last_performed_tz:
            return self.now - last_performed_tz
        else:
            return None

    @property
    def age_in_days(self):
        if age := self.age:
            return age.days
        return HealthCheckAge._NEVER_RUN_AGE

    def __str__(self):
        if not self.last_performed_tz:
            return f":dizzy_face: Never Run : {self.name}"
        emoji = HealthCheckAge._face_for_age(self.age_in_days, self.warning_age)
        return f"{emoji} {self.age_in_days} days old : {self.name}"

    _MULTIPLIER_TO_FACE = {
        0: ":simple_smile:",
        1: ":neutral_face:",
        2: ":cry:",
        3: ":rage:"
    }

    @staticmethod
    def _face_for_age(age_days: int, warning_age: timedelta):
        warning_age_days = warning_age.days
        # TODO handle warning age under 1 day
        if age_days <= 1:
            return ":smile:"
        return HealthCheckAge._MULTIPLIER_TO_FACE.get(int(age_days / warning_age_days), ":exploding_head:")

    @classmethod
    def sort_order(cls):
        return 4

    @classmethod
    def to_lines(cls, items: List['HealthCheckAge'], health_request: HealthCheckRequest) -> List[str]:
        lines: List[str] = []
        items = sorted(items, key=lambda hc: (hc.warning_age, hc.age_in_days, hc.name))
        for warning_age, warning_age_grouped in itertools.groupby(items, key=lambda hc: hc.warning_age):
            for current_age, current_and_warning_grouped in itertools.groupby(warning_age_grouped, key=lambda hc: hc.age_in_days):
                item_names = ", ".join(item.name for item in current_and_warning_grouped)
                if current_age == HealthCheckAge._NEVER_RUN_AGE:
                    lines.append(f":dizzy_face: Never Run : {item_names}")
                else:
                    if 1 < current_age <= warning_age.days:
                        pass  # don't notify if everything is within date range
                    else:
                        emoji = HealthCheckAge._face_for_age(current_age, warning_age)
                        lines.append(f"{emoji} {current_age} day{'s' if current_age != 1 else ''} old : {item_names}")
        return lines


@dataclass
class HealthCheckCustom(HealthCheckStat):
    text: str

    def __str__(self):
        return self.text

    @classmethod
    def sort_order(cls):
        return 5


health_check_signal = django.dispatch.Signal()


def populate_health_check(notification: NotificationBuilder, since: Optional[datetime] = None):
    # calling send_robust will respond with a list of tuples of receiver, the result of calling the receiver
    # the results could be an Exception, a HealthCheckStat, or a list of HealthCheckStats
    # flatten the results, report all the exceptions, then present all the HealthChecks

    now = localtime()
    if since is None:
        since = now - timedelta(days=1)
    health_request = HealthCheckRequest(since=since, now=now)

    results = []
    for caller, result in health_check_signal.send_robust(sender=None, health_request=health_request):
        if isinstance(result, Exception):
            notification.add_markdown(f"Exception generating health check by {caller}: {result}")
        else:
            results.append(result)
    checks: List[HealthCheckStat] = flatten_nested_lists(results)

    checks = sorted(checks, key=lambda hc: hc.sort_order())
    grouped_checks = [(key, list(values)) for key, values in itertools.groupby(checks, type)]
    grouped_checks = sorted(grouped_checks, key=lambda gc: gc[0].sort_order())
    recent_lines = []
    overall_lines = []
    for check_type, checks_typed in grouped_checks:
        section_lines = check_type.to_lines(checks_typed, health_request=health_request)
        if check_type.is_recent_activity():
            recent_lines.extend(section_lines)
        else:
            overall_lines.extend(section_lines)

    if recent_lines:
        notification.add_markdown("\n".join(recent_lines), indented=True)

    if overall_lines:
        notification.add_markdown("*Overall*")
        notification.add_markdown("\n".join(overall_lines), indented=True)
