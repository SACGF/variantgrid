from collections import defaultdict
from typing import Iterator

from dateutil.tz import gettz
from django.conf import settings
from django.contrib.auth.models import User

from eventlog.models import ViewEvent
from dataclasses import dataclass
from datetime import timedelta, datetime, date
from library.utils import ExportRow, export_column


@dataclass
class ReportDataRow(ExportRow):
    date: date
    user: int
    search_counts: int

    @export_column(label="Date")
    def date_column(self) -> date:
        return self.date

    @export_column(label="User")
    def user_column(self) -> int:
        return self.user

    @export_column(label="Number of Searches")
    def searches_column(self) -> int:
        return self.search_counts


def week_start_date(freq: datetime.date) -> datetime.date:
    return freq - timedelta(days=date.weekday(freq))


def stream_report_rows(interval) -> Iterator[ReportDataRow]:

    end_of_time_period = datetime.now() - timedelta(days=180)
    end_of_time_period = end_of_time_period.date()
    report_data = defaultdict(lambda: {'users': set(), 'count': 0})
    admins = User.objects.filter(is_superuser=True)
    testers_and_bots = User.objects.filter(groups__name__in={'variantgrid/tester', 'variantgrid/bot'})

    search_metric_data = ViewEvent.objects.filter(
            created__gte=end_of_time_period,
            view_name='variantopedia:search',
    ).exclude(args__search="").exclude(user__in=admins).exclude(user__in=testers_and_bots)

    ve: ViewEvent
    for ve in search_metric_data.iterator():
        created_dt: datetime = ve.created
        the_date = created_dt.astimezone(gettz(settings.TIME_ZONE)).date()
        if interval != 1:
            the_date = week_start_date(the_date)

        if user_id := ve.user_id:
            report_dict = report_data[the_date]
            report_dict["users"].add(user_id)
            report_dict["count"] += 1

    for report_date, data in report_data.items():
        user_list = len(data['users'])
        searches = data['count']
        yield ReportDataRow(date=report_date, user=user_list, search_counts=searches)


def download_search_data(request):
    frequency = 1
    if frequency_param := request.GET.get('frequency'):
        frequency = int(frequency_param)
    return ReportDataRow.streaming_csv(stream_report_rows(frequency),
                                       filename="search_data.csv")
