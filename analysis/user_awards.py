"""
Analysis award definitions (#1819): tagging and analysis work - registered from AnalysisConfig.ready()
"""
from collections import defaultdict
from datetime import datetime
from typing import Optional
from zoneinfo import ZoneInfo

from auditlog.models import LogEntry
from django.conf import settings
from django.contrib.contenttypes.models import ContentType
from django.db.models import Count
from django.db.models.functions import ExtractHour

from analysis.models import Analysis, VariantTag
from snpdb.models import UserSettingsOverride
from snpdb.models.models_enums import AwardPeriod, UserAwardKind
from snpdb.user_awards import AwardCounts, AwardDefinition, counts_by_user, register_award

ALL_PERIODS = (AwardPeriod.ALL_TIME, AwardPeriod.MONTH, AwardPeriod.DAY)


def _tags_per_tag(since: Optional[datetime]) -> AwardCounts:
    qs = VariantTag.objects.all()
    if since:
        qs = qs.filter(created__gte=since)
    return counts_by_user(qs.values("user_id", "tag_id").annotate(count=Count("id")), subject_key="tag_id")


def _tags_created(_since: Optional[datetime]) -> AwardCounts:
    return counts_by_user(VariantTag.objects.values("user_id").annotate(count=Count("id")))


def _analyses_touched(since: Optional[datetime]) -> AwardCounts:
    """ Distinct analyses with an audit log entry by the user - the Analysis itself, or one of its
        nodes (NodeAuditLogMixin puts analysis_id in additional_data) """
    analysis_ct_id = ContentType.objects.get_for_model(Analysis).pk
    qs = LogEntry.objects.filter(actor__isnull=False, content_type__app_label="analysis")
    if since:
        qs = qs.filter(timestamp__gte=since)
    per_user: dict[int, set[str]] = defaultdict(set)
    for actor_id, ct_id, object_pk, additional_data in qs.values_list("actor_id", "content_type_id", "object_pk", "additional_data").iterator():
        if ct_id == analysis_ct_id:
            analysis_id = object_pk
        elif additional_data and (analysis_id := additional_data.get("analysis_id")):
            analysis_id = str(analysis_id)
        else:
            continue
        per_user[actor_id].add(analysis_id)
    return {None: {user_id: len(analyses) for user_id, analyses in per_user.items()}}


def _tags_in_hours(hours: set[int]):
    """ Tags created during these hours in each user's own timezone """
    def counter(_since: Optional[datetime]) -> AwardCounts:
        user_timezones = dict(UserSettingsOverride.objects.exclude(timezone__isnull=True).exclude(timezone="").values_list("user_id", "timezone"))
        counts = {}
        for user_id in VariantTag.objects.values_list("user_id", flat=True).distinct():
            tz = ZoneInfo(user_timezones.get(user_id) or settings.TIME_ZONE)
            qs = VariantTag.objects.filter(user_id=user_id).annotate(hour=ExtractHour("created", tzinfo=tz))
            if count := qs.filter(hour__in=hours).count():
                counts[user_id] = count
        return {None: counts}
    return counter


register_award(AwardDefinition(
    key="top_tagger",
    kind=UserAwardKind.TITLE,
    title="Top tagger",
    description="Most variant tags created, overall and per tag",
    icon="fa-tag",
    counter=_tags_per_tag,
    periods=ALL_PERIODS,
))

register_award(AwardDefinition(
    key="top_analyst",
    kind=UserAwardKind.TITLE,
    title="Top analyst",
    description="Most analyses worked on",
    icon="fa-project-diagram",
    counter=_analyses_touched,
    periods=ALL_PERIODS,
))

register_award(AwardDefinition(
    key="tagger",
    kind=UserAwardKind.BADGE,
    title="Tagger",
    description="Variant tags created",
    icon="fa-tags",
    counter=_tags_created,
    tiers=(100, 1000, 10000),
))

register_award(AwardDefinition(
    key="analyst",
    kind=UserAwardKind.BADGE,
    title="Analyst",
    description="Analyses worked on",
    icon="fa-project-diagram",
    counter=_analyses_touched,
    tiers=(10, 100, 1000),
))

register_award(AwardDefinition(
    key="night_owl",
    kind=UserAwardKind.BADGE,
    title="Night owl",
    description="Variant tags created between 10pm and 5am",
    icon="fa-moon",
    counter=_tags_in_hours({22, 23, 0, 1, 2, 3, 4}),
    tiers=(10, 100, 1000),
    hidden=True,
))

register_award(AwardDefinition(
    key="early_bird",
    kind=UserAwardKind.BADGE,
    title="Early bird",
    description="Variant tags created between 5am and 7am",
    icon="fa-sun",
    counter=_tags_in_hours({5, 6}),
    tiers=(10, 100, 1000),
    hidden=True,
))
