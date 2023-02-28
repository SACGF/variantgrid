from typing import List

from django.dispatch import receiver

from classification.enums import ShareLevel
from classification.models import Classification, classification_flag_types, ImportedAlleleInfo, \
    ImportedAlleleInfoStatus
from flags.models import FlagType
from flags.models.flag_health_check import flag_chanced_since
from library.health_check import health_check_signal, \
    HealthCheckRequest, HealthCheckTotalAmount, HealthCheckRecentActivity, HealthCheckStat
from library.utils import pretty_label


@receiver(signal=health_check_signal)
def allele_info_health_check(sender, health_request: HealthCheckRequest, **kwargs):
    output: List[HealthCheckStat] = []

    for status in [ImportedAlleleInfoStatus.PROCESSING, ImportedAlleleInfoStatus.MATCHED_IMPORTED_BUILD]:
        not_complete = ImportedAlleleInfo.objects.filter(status=status).count()
        if not_complete:
            output.append(HealthCheckTotalAmount(
                emoji=":hourglass_flowing_sand:",
                amount=not_complete,
                name=f"Imported Allele Matching in status of \"{status.label}\""
            ))

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
            extra = ", ".join(f"{failure.imported_c_hgvs} {failure.imported_genome_build_patch_version}" for failure in last_failures)

        output.append(HealthCheckRecentActivity(
            emoji=":warning:",
            name="Variant Matching Issues",
            amount=last_failures_count,
            extra=extra,
            stand_alone=True  # doesn't make sense ot merge this with other activity numbers
        ))

    return output


@receiver(signal=health_check_signal)
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
    classifications_of_interest = Classification.dashboard_report_classifications_of_interest(since=health_request.since)
    new_classification_count = Classification.dashboard_report_new_classifications(since=health_request.since)

    return [
        HealthCheckRecentActivity(
            emoji=":blue_book:",
            amount=new_classification_count,
            name="Classifications",
            sub_type="Created"
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
