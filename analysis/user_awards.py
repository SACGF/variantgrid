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
from django.db.models import Count, Q
from django.db.models.functions import ExtractHour

from analysis.models import Analysis, VariantTag
from snpdb.models import UserSettingsOverride
from snpdb.models.models_enums import AwardPeriod, UserAwardKind
from snpdb.user_awards import AwardCounts, AwardDefinition, counts_by_user, register_award

ALL_PERIODS = (AwardPeriod.ALL_TIME, AwardPeriod.MONTH, AwardPeriod.DAY)


def _tags_created(since: Optional[datetime]) -> AwardCounts:
    qs = VariantTag.objects.all()
    if since:
        qs = qs.filter(created__gte=since)
    return counts_by_user(qs.values("user_id").annotate(count=Count("id")))


def _analyses_worked_on(since: Optional[datetime]) -> AwardCounts:
    """ Distinct analyses the user created, tagged variants in, or has an audit log entry for (the
        Analysis itself, or one of its nodes - NodeAuditLogMixin puts analysis_id in additional_data).
        Analysis auditing only started in VG4, so creation and tagging cover the history before it """
    analysis_ct_id = ContentType.objects.get_for_model(Analysis).pk
    log_qs = LogEntry.objects.filter(actor__isnull=False, content_type__app_label="analysis")
    if since:
        log_qs = log_qs.filter(timestamp__gte=since)
    logged_analysis_ids: dict[int, set[int]] = defaultdict(set)
    for actor_id, ct_id, object_pk, additional_data in log_qs.values_list("actor_id", "content_type_id", "object_pk", "additional_data").iterator():
        if ct_id == analysis_ct_id:
            logged_analysis_ids[actor_id].add(int(object_pk))
        elif additional_data and (analysis_id := additional_data.get("analysis_id")):
            logged_analysis_ids[actor_id].add(int(analysis_id))

    counts = {}
    for user_id in _users_with_analysis_activity(since, logged_analysis_ids):
        created = Q(user_id=user_id, template_type__isnull=True)
        tagged = Q(varianttag__user_id=user_id)
        if since:
            created &= Q(created__gte=since)
            tagged &= Q(varianttag__created__gte=since)
        worked_on = created | tagged | Q(pk__in=logged_analysis_ids.get(user_id, ()))
        if count := Analysis.objects.filter(worked_on).distinct().count():
            counts[user_id] = count
    return {None: counts}


def _users_with_analysis_activity(since: Optional[datetime], logged_analysis_ids: dict[int, set[int]]) -> set[int]:
    created_qs = Analysis.objects.filter(template_type__isnull=True)
    tags_qs = VariantTag.objects.filter(analysis__isnull=False)
    if since:
        created_qs = created_qs.filter(created__gte=since)
        tags_qs = tags_qs.filter(created__gte=since)
    user_ids = set(logged_analysis_ids)
    user_ids.update(created_qs.values_list("user_id", flat=True).distinct())
    user_ids.update(tags_qs.values_list("user_id", flat=True).distinct())
    return user_ids


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
    description="Most variant tags created",
    icon="fa-tag",
    counter=_tags_created,
    periods=ALL_PERIODS,
))

register_award(AwardDefinition(
    key="top_analyst",
    kind=UserAwardKind.TITLE,
    title="Top analyst",
    description="Most analyses created, tagged in or edited",
    icon="fa-project-diagram",
    counter=_analyses_worked_on,
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
    description="Analyses created, tagged in or edited",
    icon="fa-project-diagram",
    counter=_analyses_worked_on,
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
