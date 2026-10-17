"""
Classification counts per Lab: LabClassificationStats(lab) is what the lab card shows (snpdb/templatetags/lab_tags.py),
and labs_with_classifications_qs picks the labs worth listing on the labs and grouping export pages.
"Classifications" means shared ones or all, as settings.CLASSIFICATION_STATS_USE_SHARED says.
"""
from datetime import datetime
from functools import cached_property
from typing import Optional

from django.conf import settings
from django.db.models import Exists, OuterRef, QuerySet
from django.utils.timezone import now

from classification.models.classification import Classification
from snpdb.models.models import Lab
from snpdb.models.models_enums import ShareLevel


def _shared_classifications_qs() -> QuerySet[Classification]:
    return Classification.objects.filter(share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS).exclude(withdrawn=True)


def labs_with_classifications_qs(shared_only: Optional[bool] = None) -> QuerySet[Lab]:
    """ Labs with at least one classification - shared_only defaults to settings.CLASSIFICATION_STATS_USE_SHARED """
    if shared_only is None:
        shared_only = settings.CLASSIFICATION_STATS_USE_SHARED
    classifications_qs = _shared_classifications_qs() if shared_only else Classification.objects.all()
    return Lab.objects.filter(Exists(classifications_qs.filter(lab=OuterRef("pk"))))


class LabClassificationStats:

    def __init__(self, lab: Lab):
        self.lab = lab

    @cached_property
    def shared_classifications(self) -> QuerySet[Classification]:
        return _shared_classifications_qs().filter(lab=self.lab)

    @cached_property
    def classifications(self) -> QuerySet[Classification]:
        if settings.CLASSIFICATION_STATS_USE_SHARED:
            return self.shared_classifications
        return Classification.objects.filter(lab=self.lab)

    @cached_property
    def total_classifications(self) -> int:
        return self.classifications.count()

    @cached_property
    def total_shared_classifications(self) -> int:
        return self.shared_classifications.count()

    @cached_property
    def total_unshared_classifications(self) -> int:
        not_withdrawn = Classification.objects.filter(lab=self.lab).exclude(withdrawn=True).count()
        return not_withdrawn - self.total_shared_classifications

    @cached_property
    def last_created(self) -> Optional[datetime]:
        return self.classifications.order_by("created").values_list("created", flat=True).last()

    @cached_property
    def first_classification_ever_shared_date(self) -> Optional[datetime]:
        """ Includes since-withdrawn records """
        shared_qs = Classification.objects.filter(lab=self.lab, share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS)
        return shared_qs.order_by("created").values_list("created", flat=True).first()

    @cached_property
    def classifications_per_week(self) -> float:
        if first_shared := self.first_classification_ever_shared_date:
            days = max((now() - first_shared).days, 1)
            return self.total_classifications / days * 7
        return 0
