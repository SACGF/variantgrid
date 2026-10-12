"""
Classification award definitions (#1819) - registered from ClassificationConfig.ready()
"""
from datetime import timedelta

from django.db.models import Count, F, OuterRef, Subquery

from classification.models import Classification, ClassificationModification, ConditionTextMatch
from library.django_utils.model_utils import ArrayLength
from snpdb.user_awards import AwardCounts, AwardDefinition, counts_by_user, register_award

COLD_CASE_GAP = timedelta(days=183)


def _classifications_created() -> AwardCounts:
    return counts_by_user(Classification.objects.values("user_id").annotate(count=Count("id")))


def _cold_cases() -> AwardCounts:
    """ Classifications the user published a modification for when the previous modification was
        more than six months older """
    previous = ClassificationModification.objects.filter(
        classification=OuterRef("classification"), created__lt=OuterRef("created")
    ).order_by("-created").values("created")[:1]
    qs = ClassificationModification.objects.filter(published=True).annotate(
        previous_created=Subquery(previous)
    ).filter(created__gt=F("previous_created") + COLD_CASE_GAP)
    return counts_by_user(qs.values("user_id").annotate(count=Count("classification_id", distinct=True)))


def _condition_matches() -> AwardCounts:
    qs = ConditionTextMatch.objects.filter(last_edited_by__isnull=False).annotate(
        num_xrefs=ArrayLength("condition_xrefs")
    ).filter(num_xrefs__gt=0)
    return counts_by_user(qs.values(user_id=F("last_edited_by")).annotate(count=Count("id")))


register_award(AwardDefinition(
    key="cold_case",
    title="Cold case",
    description="Classifications revisited after more than six months untouched",
    icon="fa-user-secret",
    counter=_cold_cases,
    tiers=(1, 5, 25),
))

register_award(AwardDefinition(
    key="matchmaker",
    title="Matchmaker",
    description="Condition texts matched to ontology terms",
    icon="fa-link",
    counter=_condition_matches,
    tiers=(10, 100, 1000),
))

register_award(AwardDefinition(
    key="classifier",
    title="Classifier",
    description="Classifications created",
    icon="fa-clipboard",
    counter=_classifications_created,
    tiers=(10, 100, 1000),
))
