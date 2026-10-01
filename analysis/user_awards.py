"""
Analysis award definitions (#1819): tagging and analysis work - registered from AnalysisConfig.ready()
"""
from collections import defaultdict

from auditlog.models import LogEntry
from django.contrib.contenttypes.models import ContentType
from django.db.models import Count, Q

from analysis.models import Analysis, VariantTag
from snpdb.user_awards import AwardCounts, AwardDefinition, counts_by_user, register_award


def _tags_created() -> AwardCounts:
    return counts_by_user(VariantTag.objects.values("user_id").annotate(count=Count("id")))


def _analyses_worked_on() -> AwardCounts:
    """ Distinct analyses the user created, tagged variants in, or has an audit log entry for (the
        Analysis itself, or one of its nodes - NodeAuditLogMixin puts analysis_id in additional_data).
        Analysis auditing only started in VG4, so creation and tagging cover the history before it """
    analysis_ct_id = ContentType.objects.get_for_model(Analysis).pk
    log_qs = LogEntry.objects.filter(actor__isnull=False, content_type__app_label="analysis")
    logged_analysis_ids: dict[int, set[int]] = defaultdict(set)
    for actor_id, ct_id, object_pk, additional_data in log_qs.values_list("actor_id", "content_type_id", "object_pk", "additional_data").iterator():
        if ct_id == analysis_ct_id:
            logged_analysis_ids[actor_id].add(int(object_pk))
        elif additional_data and (analysis_id := additional_data.get("analysis_id")):
            logged_analysis_ids[actor_id].add(int(analysis_id))

    counts = {}
    for user_id in _users_with_analysis_activity(logged_analysis_ids):
        worked_on = (Q(user_id=user_id, template_type__isnull=True)
                     | Q(varianttag__user_id=user_id)
                     | Q(pk__in=logged_analysis_ids.get(user_id, ())))
        if count := Analysis.objects.filter(worked_on).distinct().count():
            counts[user_id] = count
    return counts


def _users_with_analysis_activity(logged_analysis_ids: dict[int, set[int]]) -> set[int]:
    user_ids = set(logged_analysis_ids)
    user_ids.update(Analysis.objects.filter(template_type__isnull=True).values_list("user_id", flat=True).distinct())
    user_ids.update(VariantTag.objects.filter(analysis__isnull=False).values_list("user_id", flat=True).distinct())
    return user_ids


register_award(AwardDefinition(
    key="tagger",
    title="Tagger",
    description="Variant tags created",
    icon="fa-tags",
    counter=_tags_created,
    tiers=(100, 1000, 10000),
))

register_award(AwardDefinition(
    key="analyst",
    title="Analyst",
    description="Analyses created, tagged in or edited",
    icon="fa-diagram-project",
    counter=_analyses_worked_on,
    tiers=(10, 100, 1000),
))
