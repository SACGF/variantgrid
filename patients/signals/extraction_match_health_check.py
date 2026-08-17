from collections import Counter

from django.dispatch import receiver

from library.health_check import HealthCheckTotalAmount, health_check_overall_stats_signal
from patients.models_enums import MatchStatus
from seqauto.models import SequencingSample
from snpdb.models import Sample


@receiver(signal=health_check_overall_stats_signal)
def extraction_match_health_check(sender, **kwargs):
    """ Parked extraction references are only a problem once they stop clearing on their own, so
        report Needs attention as the number and Pending as context """
    counts = Counter()
    for model in (Sample, SequencingSample):
        qs = model.objects.filter(extraction__isnull=True, extraction_reference__isnull=False)
        counts.update(qs.values_list("extraction_match_status", flat=True))

    needs_attention = counts[MatchStatus.NEEDS_ATTENTION]
    if not (needs_attention or counts[MatchStatus.PENDING]):
        return None
    return HealthCheckTotalAmount(
        emoji=":test_tube:",
        amount=needs_attention,
        name="unmatched extractions",
        extra=f"{counts[MatchStatus.PENDING]} still pending")
