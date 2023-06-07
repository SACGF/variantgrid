import time
from datetime import date, datetime, timedelta
from typing import Optional, Tuple, List

from dateutil import tz
from django.utils import timezone
from django.utils.timezone import localtime


def time_since(start: datetime) -> timedelta:
    end = time.time()
    return end - start


def local_date_string() -> str:
    """ Returns e.g. '2022-07-18' """
    return localtime(timezone.now()).strftime("%Y-%m-%d")


def calculate_age(born: datetime, died: Optional[datetime] = None) -> int:
    """ https://stackoverflow.com/a/9754466/295724 """
    age = None
    if born:
        if died is None:
            today = date.today()
            age = today.year - born.year - ((today.month, today.day) < (born.month, born.day))
        else:
            age = died.year - born.year - ((died.month, died.day) < (born.month, born.day))
    return age


def get_month_and_year(run_date) -> Tuple[int, int]:
    run_date_str = "%d" % int(run_date)
    parts = [run_date_str[i:i + 2] for i in range(0, len(run_date_str), 2)]
    return int(parts[1]), int(parts[0])


def year_month_string(y, m) -> str:
    return "%02d%02d" % (y, m)


def date_to_month_year_string(d) -> str:
    return month_year_string(d.year, d.month)


def month_year_string(y, m) -> str:
    return "%02d/%02d" % (m, int(str(y)[2:]))


def diff_month(d1, d2) -> int:
    return (d1.year - d2.year) * 12 + d1.month - d2.month


def get_months_since(start_month: int, start_year: int, month: int, year: int) -> int:
    return (year - start_year) * 12 + month - start_month


def month_range(start_month: int, start_year: int, offset_start: int, offset_end: int, fmt=year_month_string) -> List[str]:
    num_months = offset_end - offset_start
    labels = [fmt(start_year + (start_month + i) // 12, 1 + (start_month + i) % 12) for i in range(-1, num_months)]
    return labels


def http_header_date_now():
    return datetime.utcnow().strftime("%a, %d %b %Y %H:%M:%S GMT")


def parse_http_header_date(date_str: str):
    return datetime.strptime(date_str, "%a, %d %b %Y %H:%M:%S %Z").replace(tzinfo=tz.UTC)
