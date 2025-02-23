import logging
from dataclasses import field, dataclass
from datetime import datetime, timedelta
from functools import reduce

from threadlocals.threadlocals import get_request_variable, set_request_variable

from library.utils import time_since


def get_and_log_time_since(start: datetime, name='') -> timedelta:
    ts = time_since(start)
    if name:
        name += ': '
    logging.info("%sTotal time taken: %.2f seconds", name, ts)
    return ts


@dataclass
class DebugTime:
    description: str
    durations: list[timedelta] = field(default_factory=list)
    occurrences: int = 0

    def tick(self, duration: timedelta):
        self.durations.append(duration)
        self.occurrences += 1

    @property
    def duration(self):
        return reduce(lambda x, total: x + total, self.durations, timedelta())

    def __str__(self):
        if self.occurrences == 1:
            return f"{self.duration} - {self.description}"
        else:
            return f"{self.duration / self.occurrences} (x {self.occurrences}) - {self.description} "


class DebugTimer:

    def __init__(self):
        self.start = datetime.now()
        self.times: dict[str, DebugTime] = {}

    def tick(self, description: str):
        now = datetime.now()
        duration = now - self.start

        debug_time: DebugTime
        if existing := self.times.get(description):
            debug_time = existing
        else:
            debug_time = DebugTime(description)
            self.times[description] = debug_time

        debug_time.tick(duration)
        self.start = now

    def __str__(self):
        return "\n".join((str(debug_time) for debug_time in self.times.values()))


class NullTimer(DebugTimer):

    def tick(self, description: str):
        pass

    def __str__(self):
        return "NullTimer"


DebugTimer.NullTimer = NullTimer()


def reset_timer() -> DebugTimer:
    timer = DebugTimer()
    set_request_variable("timer", timer, use_threadlocal_if_no_request=True)
    return timer


def get_timer() -> DebugTimer:
    if timer := get_request_variable("timer", use_threadlocal_if_no_request=True):
        return timer
    # warning, asked for timer without reset timer being called
    # return null timer
    print(f"Using NullTimer - call reset_timer() if you want a timer")

    return NullTimer()