"""
Populating a classification the web created bare, so the request that made it could return straight away.
Entry point: populate_new_classification_task (@see classification.views.views.create_classification_object)
"""
from typing import Optional

import celery

from classification.autopopulate_evidence_keys.autopopulate_evidence_keys import (
    classification_complete_web_create,
    classification_populate_from_variant,
)
from classification.models import Classification, ClassificationModification
from snpdb.models import GenomeBuild


@celery.shared_task(queue='db_workers')
def populate_new_classification_task(classification_id: int, genome_build_name: str,
                                     refseq_transcript_accession: Optional[str],
                                     ensembl_transcript_accession: Optional[str],
                                     evidence: Optional[dict],
                                     copy_from_id: Optional[int],
                                     copy_gene_from_id: Optional[int]):
    """ Copy sources were permission checked by the request that queued this """
    classification = Classification.objects.get(pk=classification_id)
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    classification_populate_from_variant(classification, genome_build,
                                         refseq_transcript_accession=refseq_transcript_accession,
                                         ensembl_transcript_accession=ensembl_transcript_accession)

    copy_from = ClassificationModification.objects.get(pk=copy_from_id) if copy_from_id else None
    copy_gene_from = ClassificationModification.objects.get(pk=copy_gene_from_id) if copy_gene_from_id else None
    classification_complete_web_create(classification, classification.user, evidence=evidence,
                                       copy_from=copy_from, copy_gene_from=copy_gene_from)
