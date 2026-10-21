"""
Backfill of ImportedAlleleInfo.transcript_genome_ref_differ for records matched before the check existed (#1266).

Entry point: backfill_transcript_genome_ref_differ_task - one transcript's allele infos, queued per transcript by
the one_off_backfill_transcript_genome_ref_differ command so the transcript sequence stays cached across them.
"""
import logging

import celery

from classification.models.classification_variant_info_models import ImportedAlleleInfo
from genes.models.models_gene import TranscriptVersion, TranscriptVersionSequenceInfo


def transcript_has_stored_sequence(transcript_accession: str) -> bool:
    """ The matcher may resolve through a neighbouring version of the same transcript, so any stored version
        means the comparison is possible. With none, every row would come back unchecked - skip the lot """
    transcript_id, _ = TranscriptVersion.get_transcript_id_and_version(transcript_accession)
    return TranscriptVersionSequenceInfo.objects.filter(transcript_id=transcript_id).exists()


def backfill_transcript_genome_ref_differ(transcript_accession: str, allele_info_ids: list[int]) -> tuple[int, int]:
    """ Recomputes the flag for each allele info, re-running validation on those that now differ.
        Returns (checked, differ) counts """
    if not transcript_has_stored_sequence(transcript_accession):
        return 0, 0

    checked = differ = 0
    allele_info_qs = ImportedAlleleInfo.objects.filter(pk__in=allele_info_ids).order_by("pk")
    for allele_info in allele_info_qs.iterator():
        try:
            ref_differ = allele_info.calculate_variant_coordinate().transcript_genome_ref_differ
        except Exception:
            logging.exception("ImportedAlleleInfo %d: could not check transcript / genome reference", allele_info.pk)
            continue
        if ref_differ is None:
            continue
        checked += 1
        allele_info.transcript_genome_ref_differ = ref_differ
        if ref_differ:
            differ += 1
            allele_info.apply_validation(force_update=True)
        allele_info.save()
    return checked, differ


@celery.shared_task(queue="db_workers")
def backfill_transcript_genome_ref_differ_task(transcript_accession: str, allele_info_ids: list[int]):
    checked, differ = backfill_transcript_genome_ref_differ(transcript_accession, allele_info_ids)
    logging.info("%s: checked %d allele infos, %d differ", transcript_accession, checked, differ)
