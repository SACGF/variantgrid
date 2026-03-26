from datetime import date, datetime, timedelta, timezone
from typing import Optional

from dateutil import tz
from django.utils import timezone as django_timezone
from django.utils.timezone import localtime


def time_since(start: datetime) -> timedelta:
    return datetime.now() - start


def local_date_string() -> str:
    """ Returns e.g. '2022-07-18' """
    return localtime(django_timezone.now()).strftime("%Y-%m-%d")


def local_date_str_no_dash() -> str:
    """ Returns e.g. '20220718' """
    return localtime(django_timezone.now()).strftime("%Y%m%d")


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


def parse_yymm(run_date) -> tuple[int, int]:
    """ Parse a 4-digit YYMM value (e.g. 2201 for January 2022) into (month, year). """
    run_date_str = f"{int(run_date):d}"
    if len(run_date_str) != 4:
        raise ValueError(f"Expected 4-digit YYMM, got {run_date_str!r}")
    parts = [run_date_str[i:i + 2] for i in range(0, len(run_date_str), 2)]
    return int(parts[1]), int(parts[0])


def year_month_string(y, m) -> str:
    return f"{y:02d}{m:02d}"


def date_to_month_year_string(d) -> str:
    return month_year_string(d.year, d.month)


def month_year_string(y, m) -> str:
    return f"{m:02d}/{int(str(y)[2:]):02d}"


def diff_month(d1, d2) -> int:
    return (d1.year - d2.year) * 12 + d1.month - d2.month


def get_months_since(start_month: int, start_year: int, month: int, year: int) -> int:
    return (year - start_year) * 12 + month - start_month


def month_range(start_month: int, start_year: int, offset_start: int, offset_end: int, fmt=year_month_string) -> list[str]:
    num_months = offset_end - offset_start
    labels = [fmt(start_year + (start_month + i) // 12, 1 + (start_month + i) % 12) for i in range(-1, num_months)]
    return labels


def http_header_date_now():
    return datetime.utcnow().strftime("%a, %d %b %Y %H:%M:%S GMT")


def parse_http_header_date(date_str: str):
    return datetime.strptime(date_str, "%a, %d %b %Y %H:%M:%S %Z").replace(tzinfo=tz.UTC)


def utc_from_timestamp(ts) -> datetime:
    return datetime.fromtimestamp(ts, tz=timezone.utc)
