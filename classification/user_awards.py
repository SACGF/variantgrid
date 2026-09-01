"""
Classification award definitions (#1819) - registered from ClassificationConfig.ready()
"""
from datetime import datetime, timedelta
from typing import Optional

from django.db.models import Count, F, OuterRef, Subquery

from classification.models import Classification, ClassificationModification, ConditionTextMatch
from library.utils import ArrayLength
from snpdb.models.models_enums import AwardPeriod, UserAwardKind
from snpdb.user_awards import AwardCounts, AwardDefinition, counts_by_user, register_award

COLD_CASE_GAP = timedelta(days=183)


def _classifications_published(since: Optional[datetime]) -> AwardCounts:
    qs = ClassificationModification.objects.filter(published=True)
    if since:
        qs = qs.filter(created__gte=since)
    return counts_by_user(qs.values("user_id").annotate(count=Count("classification_id", distinct=True)))


def _classifications_created(_since: Optional[datetime]) -> AwardCounts:
    return counts_by_user(Classification.objects.values("user_id").annotate(count=Count("id")))


def _cold_cases(_since: Optional[datetime]) -> AwardCounts:
    """ Classifications the user published a modification for when the previous modification was
        more than six months older """
    previous = ClassificationModification.objects.filter(
        classification=OuterRef("classification"), created__lt=OuterRef("created")
    ).order_by("-created").values("created")[:1]
    qs = ClassificationModification.objects.filter(published=True).annotate(
        previous_created=Subquery(previous)
    ).filter(created__gt=F("previous_created") + COLD_CASE_GAP)
    return counts_by_user(qs.values("user_id").annotate(count=Count("classification_id", distinct=True)))


def _condition_matches(_since: Optional[datetime]) -> AwardCounts:
    qs = ConditionTextMatch.objects.filter(last_edited_by__isnull=False).annotate(
        num_xrefs=ArrayLength("condition_xrefs")
    ).filter(num_xrefs__gt=0)
    return counts_by_user(qs.values(user_id=F("last_edited_by")).annotate(count=Count("id")))


register_award(AwardDefinition(
    key="top_classifier",
    kind=UserAwardKind.TITLE,
    title="Top classifier",
    description="Most classifications published",
    icon="fa-clipboard",
    counter=_classifications_published,
    periods=(AwardPeriod.ALL_TIME, AwardPeriod.MONTH, AwardPeriod.DAY),
))

register_award(AwardDefinition(
    key="cold_case",
    kind=UserAwardKind.BADGE,
    title="Cold case",
    description="Classifications revisited after more than six months untouched",
    icon="fa-user-secret",
    counter=_cold_cases,
    tiers=(1, 5, 25),
))

register_award(AwardDefinition(
    key="matchmaker",
    kind=UserAwardKind.BADGE,
    title="Matchmaker",
    description="Condition texts matched to ontology terms",
    icon="fa-link",
    counter=_condition_matches,
    tiers=(10, 100, 1000),
))

register_award(AwardDefinition(
    key="classifier",
    kind=UserAwardKind.BADGE,
    title="Classifier",
    description="Classifications created",
    icon="fa-clipboard",
    counter=_classifications_created,
    tiers=(10, 100, 1000),
))
