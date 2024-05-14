from collections import defaultdict
from dataclasses import dataclass, field
from datetime import timedelta, datetime, date
from typing import Iterator
from django.contrib.auth.models import User
from django.db.models import Q
from eventlog.models import ViewEvent
from library.django_utils import require_superuser
from library.utils import ExportRow, export_column


class ReportKeys:
    SEARCHES = "searches"
    ALLELES = "alleles"
    GENE_SYMBOLS = "gene_symbols"
    CLASSIFICATIONS = "classifications"


@dataclass
class UserActivitiesRow(ExportRow):
    date: date
    unique_users: set[int] = field(default_factory=set)
    view_counts: dict = field(default_factory=lambda: defaultdict(int))
    overall_activity: int = 0

    @export_column(label="Date")
    def date_column(self) -> date:
        return self.date

    @export_column(label="Users")
    def user_column(self) -> int:
        return len(self.unique_users)

    @export_column(label="Searches")
    def searches_column(self) -> int:
        return self.view_counts[ReportKeys.SEARCHES]

    @export_column(label="Alleles Viewed")
    def alleles_viewed_column(self) -> int:
        return self.view_counts[ReportKeys.ALLELES]

    @export_column(label="Gene Symbols Viewed")
    def gene_symbols_viewed_column(self) -> int:
        return self.view_counts[ReportKeys.GENE_SYMBOLS]

    @export_column(label="Classifications Viewed")
    def classifications_viewed_column(self) -> int:
        return self.view_counts[ReportKeys.CLASSIFICATIONS]

    @export_column(label="Overall Activity")
    def overall_activity_column(self) -> int:
        return self.overall_activity


def week_start_date(freq: datetime.date) -> datetime.date:
    return freq - timedelta(days=date.weekday(freq))


_VIEW_NAMES_TO_KEYS = {
    'variantopedia:search': ReportKeys.SEARCHES,
    'variantopedia:view_allele': ReportKeys.ALLELES,
    'genes:view_gene_symbol': ReportKeys.GENE_SYMBOLS,
    'classification:view_classification': ReportKeys.CLASSIFICATIONS,
}


def stream_user_activity_rows(interval: timedelta) -> Iterator[UserActivitiesRow]:
    """
    Returns a generator of UserActivitiesRow used for making a CSV report on how much the product was used per interval
    Excludes admins, bots, testers from metrics
    :param interval: A timedelta describing how big the rows should be - should typically be 1 or 7
    """

    start_of_time_period = ViewEvent.objects.order_by('created').first().created.date()
    if interval.days > 1:
        # if we're reporting weekly, start on a Monday
        start_of_time_period = start_of_time_period - timedelta(days=start_of_time_period.weekday())
    report_end = start_of_time_period + interval
    report_data = UserActivitiesRow(
        date=start_of_time_period
    )

    # don't include tests, bots or admins in report
    excluded_users = User.objects.filter(
        Q(is_superuser=True) | Q(groups__name__in={'variantgrid/tester', 'variantgrid/bot'})
    )

    # exclude blank searches (which is the user going to the advanced search) as well as the
    # occasionally blank gene symbol when viewing a gene symbol (caused by admins putting in bad URL data typically)
    ve_qs = ViewEvent.objects.all()\
        .exclude(user__in=excluded_users)\
        .exclude(Q(view_name='variantopedia:search') & Q(args__search=""))\
        .exclude(Q(view_name='genes:view_gene_symbol') & Q(args__gene_symbol="")) \
        .order_by('created')

    ve: ViewEvent
    for ve in ve_qs.iterator():
        while ve.created.date() >= report_end:
            # cycle through report data dates until we have a report data interval for the current view event
            yield report_data
            report_data = UserActivitiesRow(date=report_end)
            report_end = report_end + interval

        # if the view event is for a specific page type, increment the count
        if view_key := _VIEW_NAMES_TO_KEYS.get(ve.view_name):
            report_data.view_counts[view_key] += 1

        # if the view had a user (nearly all of them should) record the user_id and increment overall activity
        if user_id := ve.user_id:
            report_data.overall_activity += 1
            report_data.unique_users.add(user_id)
    yield report_data


@require_superuser
def download_search_data(request):
    frequency = 1
    if frequency_param := request.GET.get('frequency'):
        frequency = int(frequency_param)
    return UserActivitiesRow.streaming_csv(
        stream_user_activity_rows(timedelta(days=frequency)),
        filename="search_data.csv"
    )
