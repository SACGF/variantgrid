"""
Independent feeds arrive in any order: the client may post a run before it accessions the specimen, or
accession it a day later. A row whose extraction cannot be resolved yet is parked rather than rejected,
and this re-resolves the parked ones - on a schedule, and again whenever new extractions land.

Entry point is reconcile_pending_extractions, over every ExtractionMatchMixin model that parks a claim
(SequencingSample, Sample) and over LibraryQC, whose claim is on the Specimen its pair was taken from
and whose arm is linked to its sample sheet row once the run's sheet is registered.
"""
from datetime import timedelta

import celery
from django.conf import settings
from django.db.models import Q
from django.utils import timezone

from patients.external_references import ExternalReference, resolve_reference
from patients.models import Extraction, Specimen
from patients.models_enums import MatchStatus
from seqauto.models import (
    SAMPLE_SHEET_PAIR_ID_COLUMN,
    LibraryQC,
    SequencingRun,
    SequencingSample,
    sequencing_sample_for_pair,
)
from snpdb.models import Sample

PENDING_STATES = [MatchStatus.PENDING, MatchStatus.NEEDS_ATTENTION]

_COUNT_KEYS = {
    MatchStatus.MATCHED: "matched",
    MatchStatus.PENDING: "still_pending",
    MatchStatus.NEEDS_ATTENTION: "needs_attention",
}


@celery.shared_task(queue='db_workers')
def reconcile_pending_extractions() -> dict:
    counts = {"matched": 0, "still_pending": 0, "needs_attention": 0, "from_sequencing_sample": 0,
              "library_qc_linked": 0}

    for model in (SequencingSample, Sample):
        qs = model.objects.filter(extraction__isnull=True,
                                  extraction_reference__isnull=False,
                                  extraction_match_status__in=PENDING_STATES)
        for row in qs.iterator():
            user = _row_user(row)
            reference = ExternalReference.from_data(row.extraction_reference)
            resolved = resolve_reference(Extraction, reference, user)
            if resolved.status == MatchStatus.PENDING and _past_pending_window(row.extraction_match_date):
                # Past the window this is a real mismatch rather than the load race, and wants a human
                resolved.status = MatchStatus.NEEDS_ATTENTION
            row.apply_extraction_match(resolved)
            counts[_COUNT_KEYS[resolved.status]] += 1

    # A MetricsOutput can land before the CombinedVariantOutput accessions the case it describes
    library_qc = LibraryQC.objects.filter(specimen__isnull=True,
                                          specimen_match_status__in=PENDING_STATES).exclude(specimen_reference="")
    for row in library_qc.iterator():
        reference = ExternalReference(reference_id=row.specimen_reference)
        resolved = resolve_reference(Specimen, reference, row.user)
        if resolved.status == MatchStatus.PENDING and _past_pending_window(row.specimen_match_date):
            resolved.status = MatchStatus.NEEDS_ATTENTION
        row.apply_specimen_match(resolved)
        counts[_COUNT_KEYS[resolved.status]] += 1

    counts["library_qc_linked"] = link_library_qc_to_sequencing_samples()

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


def link_library_qc_to_sequencing_samples() -> int:
    """ A MetricsOutput can land before its run's sample sheet, and a re-sent sheet builds new
        SequencingSample rows: rows with no arm, or an arm off a superseded sheet, are looked up again
        on every run whose current sheet carries the Pair_ID column. The run link is filled in at the
        same time where the file arrived before the run was registered """
    runs_with_pairs = SequencingRun.objects.filter(
        sequencingruncurrentsamplesheet__sample_sheet__sequencingsample__sequencingsampledata__column=SAMPLE_SHEET_PAIR_ID_COLUMN,
    ).distinct()
    stale = LibraryQC.objects.filter(sequencing_run_name__in=runs_with_pairs.values("name")).filter(
        Q(sequencing_sample__isnull=True) |
        Q(sequencing_sample__sample_sheet__sequencingruncurrentsamplesheet__isnull=True))
    linked = 0
    for row in stale.iterator():
        if row.sequencing_run_id is None:
            row.sequencing_run = SequencingRun.objects.get(name=row.sequencing_run_name)
        sequencing_sample = sequencing_sample_for_pair(row.sequencing_run, row.pair_id, row.nucleic_acid)
        if sequencing_sample is None:
            continue
        row.sequencing_sample = sequencing_sample
        row.save()
        linked += 1
    return linked


def _row_user(row):
    """ Who the claim is re-resolved as - the row's own uploader where it has one, else the owner of
        the VCF the sample came in on """
    if user := getattr(row, "user", None):
        return user
    return getattr(getattr(row, "vcf", None), "user", None)


def _past_pending_window(parked) -> bool:
    if parked is None:
        return False
    days = settings.PATIENT_EXTRACTION_MATCH_PENDING_DAYS
    return timezone.now() - parked > timedelta(days=days)
