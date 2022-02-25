from collections import defaultdict
from dataclasses import dataclass
from datetime import timedelta, datetime
from typing import Any, List, Callable, TypeVar, Generic, Optional

from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.db.models import Model, Count
from django.http import HttpRequest
from django.http.response import HttpResponseBase
from django.shortcuts import render
from django.utils.timezone import now
from lazy import lazy
from classification.models import Classification
from eventlog.models import ViewEvent
from snpdb.models import Allele

T = TypeVar("T")

@dataclass
class Counted(Generic[T]):
    pk: Any
    count: int
    resolver: Optional[Callable[[Any], T]]

    @lazy
    def resolve(self):
        if resolver := self.resolver:
            return resolver(self.pk)
        return self.pk

    def __lt__(self, other):
        return self.count < other.count


@dataclass(frozen=True)
class ViewEventCounts:
    time_ago: timedelta

    @lazy
    def as_of(self) -> datetime:
        return now() - self.time_ago

    def count_field(self, field_name: str, resolver: Optional[Callable]) -> List[Counted]:
        id_to_count = defaultdict(int)
        field_name = f"args__{field_name}"
        for pk, user in ViewEvent.objects.filter(**{f"{field_name}__isnull": False}).filter(created__gte=self.as_of).annotate(
                total=Count(field_name)).order_by(field_name).values_list(field_name, 'user'):
            id_to_count[pk] += 1  # just want to count once per unique user

        return sorted((Counted(pk, count, resolver) for pk, count in id_to_count.items()), reverse=True)

    @staticmethod
    def resolver_for_model(model: Model):
        def resolver(pk: Any):
            return model.objects.filter(pk=pk).first()
        return resolver

    @lazy
    def classification_views(self) -> List[Counted[Classification]]:
        return self.count_field("classification_id", ViewEventCounts.resolver_for_model(Classification))

    @lazy
    def allele_views(self) -> List[Counted[Allele]]:
        return self.count_field("allele_id", ViewEventCounts.resolver_for_model(Allele))

    @lazy
    def page_views(self) -> List[Counted[str]]:
        id_to_count = defaultdict(int)
        for pk in ViewEvent.objects.filter(
                created__gte=self.as_of).annotate(
                total=Count('user')).order_by('user').values_list('view_name', flat=True):
            id_to_count[pk] += 1  # just want to count once per unique user

        return sorted((Counted(pk, count, None) for pk, count in id_to_count.items()), reverse=True)

    @lazy
    def active_users(self) -> List[Counted[str]]:
        id_to_count = defaultdict(int)
        for pk, total in ViewEvent.objects.filter(
                created__gte=self.as_of).annotate(
            total=Count('user')).order_by('user').values_list('user', 'total'):
            id_to_count[pk] += total  # just want to count once per unique user

        return sorted((Counted(pk, count, ViewEventCounts.resolver_for_model(User)) for pk, count in id_to_count.items()), reverse=True)

def view_classifiaction_metrics(request: HttpRequest) -> HttpResponseBase:
    if not request.user.is_superuser:
        raise PermissionDenied()

    context = {
        "counts": ViewEventCounts(time_ago=timedelta(30))
    }
    return render(request, "classification/classification_view_metrics.html", context)

