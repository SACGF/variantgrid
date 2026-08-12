"""
Independent feeds arrive in any order: the client may post a run before it accessions the specimen, or
accession it a day later. A row whose extraction cannot be resolved yet is parked rather than rejected,
and this re-resolves the parked ones - on a schedule, and again whenever new extractions land.
"""
from datetime import timedelta

import celery
from django.conf import settings
from django.utils import timezone

from patients.external_references import ExternalReference, resolve_reference
from patients.models import Extraction
from patients.models_enums import MatchStatus
from seqauto.models import SequencingSample
from snpdb.models import Sample

PENDING_STATES = [MatchStatus.PENDING, MatchStatus.NEEDS_ATTENTION]

_COUNT_KEYS = {
    MatchStatus.MATCHED: "matched",
    MatchStatus.PENDING: "still_pending",
    MatchStatus.NEEDS_ATTENTION: "needs_attention",
}


@celery.shared_task(queue='db_workers')
def reconcile_pending_extractions() -> dict:
    counts = {"matched": 0, "still_pending": 0, "needs_attention": 0, "from_sequencing_sample": 0}

    for model in (SequencingSample, Sample):
        qs = model.objects.filter(extraction__isnull=True,
                                  extraction_reference__isnull=False,
                                  extraction_match_status__in=PENDING_STATES)
        for row in qs.iterator():
            user = getattr(getattr(row, "vcf", None), "user", None)
            reference = ExternalReference.from_data(row.extraction_reference)
            resolved = resolve_reference(Extraction, reference, user)
            if resolved.status == MatchStatus.PENDING and _past_pending_window(row):
                # Past the window this is a real mismatch rather than the load race, and wants a human
                resolved.status = MatchStatus.NEEDS_ATTENTION
            row.apply_extraction_match(resolved)
            counts[_COUNT_KEYS[resolved.status]] += 1

    # Route 1 arriving after the VCF: the link call set SequencingSample.extraction, but
    # link_samples_and_vcfs_to_sequencing had already run and had nothing to carry down
    unlinked = Sample.objects.filter(extraction__isnull=True,
                                     samplefromsequencingsample__sequencing_sample__extraction__isnull=False)
    for sample in unlinked.select_related("samplefromsequencingsample__sequencing_sample").iterator():
        sample.extraction = sample.samplefromsequencingsample.sequencing_sample.extraction
        sample.extraction_match_status = MatchStatus.MATCHED
        sample.extraction_match_date = timezone.now()
        sample.save()
        counts["from_sequencing_sample"] += 1

    return counts


def _past_pending_window(row) -> bool:
    parked = row.extraction_match_date
    if parked is None:
        return False
    days = settings.PATIENT_EXTRACTION_MATCH_PENDING_DAYS
    return timezone.now() - parked > timedelta(days=days)
