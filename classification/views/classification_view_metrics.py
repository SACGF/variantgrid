import operator
from collections import defaultdict
from dataclasses import dataclass
from datetime import timedelta, datetime
from functools import cached_property, reduce
from typing import Any, List, Callable, TypeVar, Generic, Optional

from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.db.models import Model, Count, Q, QuerySet
from django.http import HttpRequest
from django.http.response import HttpResponseBase
from django.shortcuts import render
from django.utils.timezone import now

from classification.models import Classification, DiscordanceReport
from eventlog.models import ViewEvent
from genes.models import GeneSymbol
from library.django_utils import require_superuser
from snpdb.models import Allele

T = TypeVar("T")

@dataclass
class Counted(Generic[T]):
    pk: Any
    count: int
    resolver: Optional[Callable[[Any], T]]

    @cached_property
    def resolve(self):
        if resolver := self.resolver:
            return resolver(self.pk)
        return self.pk

    def __lt__(self, other):
        return self.count < other.count


@dataclass(frozen=True)
class ViewMetricsType:
    counts: List[Counted[Any]]
    name: str
    model_id: str


@dataclass(frozen=True)
class ViewEventCounts:
    time_ago: timedelta
    exclude_admin: bool = False

    @property
    def page_suffix(self):
        parts = [f"over the last {self.time_ago.days} days"]
        if self.exclude_admin:
            parts.append("excluding admins")
        return ", ".join(parts)

    @staticmethod
    def from_request(request: HttpRequest):
        time_ago = timedelta(days=30)
        if days_old_str := request.GET.get('days'):
            time_ago = timedelta(days=int(days_old_str))

        # default to exclude admin
        exclude_admin = request.GET.get('exclude_admin') == "true" or not request.GET.get('exclude_admin')
        return ViewEventCounts(time_ago=time_ago, exclude_admin=exclude_admin)

    @cached_property
    def as_of(self) -> datetime:
        return now() - self.time_ago

    @property
    def base_filter(self) -> Q:
        qs: List[Q] = []
        if self.exclude_admin:
            qs.append(Q(user__is_superuser=False))
        qs.append(Q(created__gte=self.as_of))
        return reduce(operator.and_, qs)

    def count_field(self, field_name: str, resolver: Optional[Callable]) -> List[Counted]:
        id_to_count = defaultdict(int)
        field_name = f"args__{field_name}"
        for values in ViewEvent.objects \
                .values(field_name, 'user')\
                .filter(self.base_filter)\
                .filter(**{f"{field_name}__isnull": False})\
                .annotate(total=Count('pk'))\
                .order_by():
            if non_blank := values.get(field_name):
                id_to_count[non_blank] += 1  # just want to count once per unique user
            # TODO limit the number of results we look at? (though wont work if we want unique users)

        return sorted((Counted(pk, count, resolver) for pk, count in id_to_count.items()), reverse=True)

    @staticmethod
    def resolver_for_model(model: Model):
        def resolver(pk: Any):
            if pk:
                if first := model.objects.filter(pk=pk).first():
                    return first
                else:
                    return f"{pk}"
            else:
                return "blank"
        return resolver

    @cached_property
    def classification_views(self) -> List[Counted[Classification]]:
        return self.count_field("classification_id", ViewEventCounts.resolver_for_model(Classification))

    @cached_property
    def discordance_report_views(self) -> List[Counted[Classification]]:
        return self.count_field("discordance_report_id", ViewEventCounts.resolver_for_model(DiscordanceReport))

    @cached_property
    def allele_views(self) -> List[Counted[Allele]]:
        return self.count_field("allele_id", ViewEventCounts.resolver_for_model(Allele))

    @cached_property
    def gene_symbol_views(self) -> List[Counted[GeneSymbol]]:
        return self.count_field("gene_symbol", ViewEventCounts.resolver_for_model(GeneSymbol))

    @cached_property
    def page_views(self) -> List[Counted[str]]:
        id_to_count = defaultdict(int)
        for values in ViewEvent.objects\
                .values('view_name', 'user')\
                .filter(self.base_filter)\
                .annotate(total=Count('pk'))\
                .order_by():
            id_to_count[values.get('view_name')] += 1  # just want to count once per unique user

        return sorted((Counted(pk, count, None) for pk, count in id_to_count.items()), reverse=True)

    @cached_property
    def active_users(self) -> List[Counted[str]]:
        id_to_count = defaultdict(int)
        for values in ViewEvent.objects \
                .values('user')\
                .filter(self.base_filter)\
                .annotate(total=Count('pk'))\
                .order_by():
            id_to_count[values.get('user')] += values.get('total')  # just want to count once per unique user

        return sorted((Counted(pk, count, ViewEventCounts.resolver_for_model(User)) for pk, count in id_to_count.items()), reverse=True)

    def view_metrics(self) -> List[ViewMetricsType]:
        return [
            ViewMetricsType(
                counts=self.gene_symbol_views,
                name="Gene Symbol",
                model_id="gene_symbol"
            ),
            ViewMetricsType(
                counts=self.allele_views,
                name="Alleles",
                model_id="allele_id"
            ),
            ViewMetricsType(
                counts=self.classification_views,
                name="Classifications",
                model_id="classification_id"
            ),
            ViewMetricsType(
                counts=self.discordance_report_views,
                name="Discordance Reports",
                model_id="discordance_report_id"
            )
        ]

    def recent_views(self) -> QuerySet[ViewEvent]:
        return ViewEvent.objects.filter(self.base_filter).order_by('-created')

    def all_views_for(self, request: HttpRequest) -> QuerySet[ViewEvent]:
        views = ViewEvent.objects \
            .filter(self.base_filter) \
            .order_by('-created') \
            .select_related('user')

        if view_name := request.GET.get('view_name'):
            views = views.filter(view_name=view_name)
        elif classification_id := request.GET.get('classification_id'):
            views = views.filter(args__classification_id=int(classification_id))
        elif allele_id := request.GET.get('allele_id'):
            views = views.filter(args__allele_id=int(allele_id))
        elif discordance_report_id := request.GET.get('discordance_report_id'):
            views = views.filter(args__discordance_report_id=int(discordance_report_id))
        elif gene_symbol_id := request.GET.get('gene_symbol'):
            views = views.filter(args__gene_symbol=gene_symbol_id)
        elif user_id := request.GET.get('user_id'):
            views = views.filter(user=user_id)

        return views


@require_superuser
def view_classification_metrics(request: HttpRequest) -> HttpResponseBase:
    if not request.user.is_superuser:
        raise PermissionDenied()

    vec = ViewEventCounts.from_request(request)

    context = {
        "counts": vec,
        "days": vec.time_ago.days,
        "exclude_admin": vec.exclude_admin,
        "days_options": [1, 7, 30, 60, 90],
        "page_suffix": vec.page_suffix
    }
    return render(request, "classification/classification_view_metrics.html", context)


@require_superuser
def view_page_metrics_detail(request: HttpRequest) -> HttpResponseBase:
    vec = ViewEventCounts.from_request(request)
    views = vec.all_views_for(request)

    context = {
        "views": views
    }
    return render(request, "classification/classification_view_metrics_detail.html", context)
