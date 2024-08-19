from dataclasses import dataclass
from enum import Enum, IntEnum
from functools import cached_property
from typing import Optional
import django
from django.contrib.auth.models import User

from library.utils import text_utils


###
# Uptime Check is similar to health check but
# The results is to be displayed to a public web page, so it needs to not show sensitive information
# The results are more about if critical systems are currently working or not, so we can get an alert from uptime robot if they're not
###

class UptimeCheckStatus(IntEnum):
    OKAY = 1
    NON_CRITICAL_FAILURE = 2
    CRITICAL_FAILURE = 3
    EXCEPTION_WHEN_CHECKING = 4

    @property
    def pretty_label(self):
        return text_utils.pretty_label(self.name)

    def severity(self) -> str:
        match self:
            case UptimeCheckStatus.OKAY: return "S"
            case UptimeCheckStatus.NON_CRITICAL_FAILURE: return "E"
            case _: return "C"

@dataclass
class UptimeCheckResponse:
    name: str
    status: UptimeCheckStatus
    note: Optional[str] = None


@dataclass
class UptimeCheckOverall:
    uptime_checks: list[UptimeCheckResponse]

    @cached_property
    def status(self) -> UptimeCheckStatus:
        return max(uc.status for uc in self.uptime_checks)


uptime_check_signal = django.dispatch.Signal()


def retrieve_uptime_response():
    all_uptimes = []

    status = UptimeCheckStatus.OKAY
    try:
        # note that it's quite likely if the database is down, something well before this will cause an exception
        # but feels dishonest to say the database is okay without checking anything
        User.objects.first()
    except:
        status = UptimeCheckStatus.CRITICAL_FAILURE
    all_uptimes.append(UptimeCheckResponse(
        name="Database",
        status=status
    ))

    for caller, result in uptime_check_signal.send_robust(sender=None):
        if result:
            if isinstance(result, Exception):
                all_uptimes.append(UptimeCheckResponse(
                    name=str(caller),
                    status=UptimeCheckStatus.EXCEPTION_WHEN_CHECKING,
                    note=f"Exception generating health check by {caller}: {result}"
                ))
            else:
                all_uptimes.append(result)

    return UptimeCheckOverall(uptime_checks=all_uptimes)
