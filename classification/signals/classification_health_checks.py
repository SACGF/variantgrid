from datetime import timedelta

from django.db.models.query_utils import Q
from django.dispatch import receiver
from django.utils import timezone as django_timezone

from classification.enums.classification_enums import ShareLevel
from classification.models import (
    Classification,
    ImportedAlleleInfo,
    ImportedAlleleInfoStatus,
    classification_flag_types,
)
from flags.models import Flag, FlagStatus, FlagType
from flags.models.flag_health_check import flag_chanced_since
from library.health_check import (
    HealthCheckRecentActivity,
    HealthCheckRequest,
    HealthCheckStat,
    HealthCheckTotalAmount,
    health_check_overall_stats_signal,
    health_check_signal,
)
from library.utils import limit_str

"""
Reports information about classifications to the Slack health report
"""


class ClassificationOutstandingIssues:

    def __init__(self, classification: 'Classification'):
        self.classification = classification
        self.flags = []
        self.issues = []

    @property
    def pk(self):
        return self.classification.pk

    def add_flag(self, flag: str):
        self.flags.append(flag)

    def add_issue(self, issue: str):
        self.issues.append(issue)

    def __str__(self):
        return f"({self.classification.friendly_label}) {', '.join(self.issues)} {', '.join(self.flags)}"


def dashboard_report_new_classifications(since) -> int:
    return Classification.objects.filter(created__gte=since).count()


def dashboard_report_classifications_of_interest(since) -> list[ClassificationOutstandingIssues]:
    min_age = django_timezone.now() - timedelta(minutes=2)  # give records 2 minutes to matching properly before reporting

    time_range_q = Q(created__gte=since) & Q(created__lte=min_age)

    # want to find new tags that are still open
    flag_collections = Flag.objects.filter(time_range_q, resolution__status=FlagStatus.OPEN)
    flag_collections = flag_collections.order_by('collection__id').values_list('collection__id', flat=True)
    flag_q = Q(flag_collection_id__in=flag_collections.distinct())
    missing_chgvs_q = (Q(allele_info__grch37__resolved_hgvs__isnull=True) | Q(allele_info__grch38__resolved_hgvs__isnull=True))
    coi_qs = Classification.objects.filter(flag_q | (time_range_q & missing_chgvs_q))
    coi_qs = coi_qs.order_by('-pk').select_related('lab', 'flag_collection')

    summaries: list[ClassificationOutstandingIssues] = []
    c: Classification
    for c in coi_qs:
        coi = ClassificationOutstandingIssues(c)
        this_flags = Flag.objects.filter(time_range_q, resolution__status=FlagStatus.OPEN,
                                         collection=c.flag_collection)

        variant_matching = False
        for flag_type in this_flags.order_by('flag_type').values_list('flag_type', flat=True):
            if flag_type == 'classification_matching_variant':
                variant_matching = True
            coi.add_flag(flag_type)
        if not variant_matching:
            if not c.resolved_hgvs_grch37:
                coi.add_issue("No cached 37 representation")
            if not c.resolved_hgvs_grch38:
                coi.add_issue("No cached 38 representation")

        summaries.append(coi)

    return summaries


@receiver(signal=health_check_signal)
def allele_info_health_check(sender, health_request: HealthCheckRequest, **kwargs):
    output: list[HealthCheckStat] = []

    last_failures = ImportedAlleleInfo.objects.filter(
        latest_validation__created__gte=health_request.since,
        latest_validation__created__lte=health_request.now,
        latest_validation__include=False
    ).\
        select_related('imported_genome_build_patch_version')
    last_failures_count = last_failures.count()

    if last_failures_count:
        extra = None
        if last_failures_count <= 5:
            extra = ", ".join(f"{limit_str(failure.imported_c_hgvs, 80)} {failure.imported_genome_build_patch_version}" for failure in last_failures)

        output.append(HealthCheckRecentActivity(
            emoji=":warning:",
            name="Variant Matching Issues",
            amount=last_failures_count,
            extra=extra,
            stand_alone=True  # doesn't make sense ot merge this with other activity numbers
        ))

    return output


@receiver(signal=health_check_overall_stats_signal)
def allele_info_overall_stats(sender, **kwargs):
    output: list[HealthCheckStat] = []

    for status in [ImportedAlleleInfoStatus.PROCESSING, ImportedAlleleInfoStatus.MATCHED_IMPORTED_BUILD]:
        not_complete = ImportedAlleleInfo.objects.filter(status=status).count()
        if not_complete:
            output.append(HealthCheckTotalAmount(
                emoji=":hourglass_flowing_sand:",
                amount=not_complete,
                name=f"Imported Allele Matching in status of \"{status.label}\""
            ))

    return output


@receiver(signal=health_check_overall_stats_signal)
def classifications_health_check_count(sender, health_request: HealthCheckRequest, **kwargs):
    total_classification_qs = Classification.objects.filter(lab__external=False, withdrawn=False, created__lte=health_request.now).exclude(lab__name__icontains='legacy')
    total = total_classification_qs.count()
    total_shared = total_classification_qs.filter(share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS).count()
    if total:
        percent_shared = 100.0 * (float(total_shared) / float(total))

        return HealthCheckTotalAmount(
            emoji=":blue_book:",
            amount=total,
            name="Classifications",
            extra=f"{int(percent_shared)}% shared"
        )


@receiver(signal=health_check_signal)
def classification_health_check_activity(sender, health_request: HealthCheckRequest, **kwargs):
    # TODO, also want to restrict this data based on now (not too much of an issue if now() is pretty much now)
    classifications_of_interest = dashboard_report_classifications_of_interest(since=health_request.since)
    new_classification_count = dashboard_report_new_classifications(since=health_request.since)

    classification_new = Classification.objects.filter(created__gte=health_request.since)

    return [
        HealthCheckRecentActivity(
            emoji=":blue_book:",
            amount=new_classification_count,
            name="Classifications",
            sub_type="Created",
            preview=[classification.preview for classification in classification_new],
        ),
        HealthCheckRecentActivity(
            emoji=":blue_book:",
            amount=len(classifications_of_interest),
            name="Classifications",
            sub_type="Of Interest"
        )
    ]


@receiver(signal=health_check_signal)
def classification_flag_health_check_activity(sender, health_request: HealthCheckRequest, **kwargs):
    flag_deltas = flag_chanced_since(
        health_request.since,
        flag_types=FlagType.objects.filter(context=classification_flag_types.classification_flag_context)
    )
    checks = []
    for flag_delta in flag_deltas:

        parts = []
        if flag_delta.added:
            parts.append(f"+{flag_delta.added}")
        if flag_delta.resolved:
            parts.append(f"-{flag_delta.resolved}")

        if parts:
            checks.append(
                HealthCheckRecentActivity(
                    emoji=":flags:",
                    name=f"Flag {flag_delta.flag_type}",
                    amount=" / ".join(parts)
                )
            )
    return checks
