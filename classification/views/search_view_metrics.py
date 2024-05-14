from collections import defaultdict
from dataclasses import dataclass
from datetime import timedelta, datetime, date
from typing import Iterator

from dateutil.tz import gettz
from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import Q, Min

from eventlog.models import ViewEvent
from library.django_utils import require_superuser
from library.utils import ExportRow, export_column


@dataclass
class ReportDataRow(ExportRow):
    date: date
    users: int
    search_counts: int
    allele_counts: int
    gene_symbols_counts: int
    classification_counts: int
    overall_activity: int

    @export_column(label="Date")
    def date_column(self) -> date:
        return self.date

    @export_column(label="Users")
    def user_column(self) -> int:
        return self.users

    @export_column(label="Number of Searches")
    def searches_column(self) -> int:
        return self.search_counts

    @export_column(label="Alleles Viewed")
    def alleles_viewed_column(self) -> int:
        return self.allele_counts

    @export_column(label="Gene Symbols Viewed")
    def gene_symbols_viewed_column(self) -> int:
        return self.gene_symbols_counts

    @export_column(label="Classifications Viewed")
    def classifications_viewed_column(self) -> int:
        return self.classification_counts

    @export_column(label="Overall Activity")
    def overall_activity_column(self) -> int:
        return self.overall_activity


def week_start_date(freq: datetime.date) -> datetime.date:
    return freq - timedelta(days=date.weekday(freq))


def stream_report_rows(interval) -> Iterator[ReportDataRow]:

    # start_of_time_period = datetime.now() - timedelta(days=365)
    latest_created_date = ViewEvent.objects.aggregate(Min('created'))['created__min']
    start_of_time_period = latest_created_date.date()
    report_data = defaultdict(lambda: {'users': set(), 'search_count': 0, 'allele_count': 0,
                                       'gene_symbol_count': 0, 'classification_count': 0,
                                       'overall_activity': 0})
    excluded_users = User.objects.filter(
        Q(is_superuser=True) | Q(groups__name__in={'variantgrid/tester', 'variantgrid/bot'}))

    metrics_definitions = [
        ('variantopedia:search', 'search_count', ~Q(args__search="")),
        ('variantopedia:view_allele', 'allele_count'),
        ('genes:view_gene_symbol', 'gene_symbol_count', ~Q(args__gene_symbol="")),
        ('classification:view_classification', 'classification_count'),
        ('', 'overall_activity'),
    ]

    ve: ViewEvent
    for view_name, key, *extra_filters in metrics_definitions:
        if view_name == '':
            metric_data = ViewEvent.objects.filter(
                created__gte=start_of_time_period
            ).exclude(user__in=excluded_users)
        else:
            metric_data = ViewEvent.objects.filter(
                created__gte=start_of_time_period,
                view_name=view_name
            ).exclude(user__in=excluded_users)

        if extra_filters:
            metric_data = metric_data.filter(*extra_filters)
        for ve in metric_data.iterator():
            created_dt: datetime = ve.created
            the_date = created_dt.astimezone(gettz(settings.TIME_ZONE)).date()
            if interval != 1:
                the_date = week_start_date(the_date)

            if user_id := ve.user_id:
                report_dict = report_data[the_date]
                report_dict["users"].add(user_id)
                report_dict[key] += 1

    current_date = week_start_date(start_of_time_period)
    today = datetime.now().date()

    while current_date <= today:
        data = report_data.get(current_date, {'users': set(), 'search_count': 0, 'allele_count': 0,
                                              'gene_symbol_count': 0, 'classification_count': 0,
                                              'overall_activity': 0})
        user_count = len(data['users'])
        searches = data['search_count']
        allele_count = data['allele_count']
        gene_symbol_count = data['gene_symbol_count']
        classification_count = data['classification_count']
        overall_activity = data['overall_activity']
        yield ReportDataRow(date=current_date, users=user_count, search_counts=searches,
                            allele_counts=allele_count, gene_symbols_counts=gene_symbol_count,
                            classification_counts=classification_count,
                            overall_activity=overall_activity)
        current_date += timedelta(days=(1 if interval == 1 else 7))


@require_superuser
def download_search_data(request):
    frequency = 1
    if frequency_param := request.GET.get('frequency'):
        frequency = int(frequency_param)
    return ReportDataRow.streaming_csv(stream_report_rows(frequency),
                                       filename="search_data.csv")
