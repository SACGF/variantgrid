import celery
from django.core.exceptions import ObjectDoesNotExist
from django.db.models.query_utils import Q
import logging
import time

from genes.canonical_transcripts.canonical_transcript_manager import CanonicalTranscriptManager
from genes.gene_matching import GeneSymbolMatcher
from genes.models import GeneCoverageCollection, GeneCoverageCanonicalTranscript
from seqauto.models import EnrichmentKit
from snpdb.models import DataState


@celery.task
def reload_gene_coverage_collection(gene_coverage_collection_id):
    logging.info("reload_gene_coverage_collection(%s) START", gene_coverage_collection_id)

    start = time.time()

    gene_coverage_collection = GeneCoverageCollection.objects.get(pk=gene_coverage_collection_id)
    gene_coverage_collection.genecoverage_set.all().delete()
    gene_coverage_collection.genecoveragecanonicaltranscript_set.all().delete()
    gene_coverage_collection.data_state = DataState.RUNNING
    gene_coverage_collection.save()

    gene_matcher = GeneSymbolMatcher()
    canonical_transcript_manager = CanonicalTranscriptManager()

    try:
        enrichment_kit = gene_coverage_collection.qcgenecoverage.qc.sequencing_sample.enrichment_kit
    except ObjectDoesNotExist:
        enrichment_kit = None

    gene_coverage_collection.load_from_file(enrichment_kit, gene_matcher=gene_matcher, canonical_transcript_manager=canonical_transcript_manager)
    gene_coverage_collection.data_state = DataState.COMPLETE
    gene_coverage_collection.save()

    end = time.time()
    logging.info("reload_gene_coverage_collection(%s) DONE in %.1f seconds", gene_coverage_collection_id, (end - start))


# TODO: This is only needed to migrate existing data - it just takes hours so want to spread across celery tasks
# Once all environments https://github.com/SACGF/variantgrid/wiki/Upgrade_Notes have this applied:
# https://github.com/SACGF/variantgrid/issues/1216#issuecomment-440561628 delete this task etc.
@celery.task
def create_canonical_gene_coverage_for_enrichment_kit(enrichment_kit_id):
    #logging.info("create_canonical_gene_coverage_for_enrichment_kit %s", enrichment_kit_id)

    canonical_transcript_manager = CanonicalTranscriptManager()

    if enrichment_kit_id:
        enrichment_kit = EnrichmentKit.objects.get(pk=enrichment_kit_id)
        canonical_collection = canonical_transcript_manager.get_canonical_collection_for_enrichment_kit(enrichment_kit)
        coverage_collection_qs = GeneCoverageCollection.objects.filter(qc__bam_file__unaligned_reads__sequencing_sample__enrichment_kit=enrichment_kit)
    else:
        canonical_collection = canonical_transcript_manager.get_default_canonical_collection()
        coverage_collection_qs = GeneCoverageCollection.objects.filter(qc__isnull=True)

    canonical_transcripts = canonical_transcript_manager.get_canonical_transcripts(canonical_collection)

    # Skip ones that have already been calculated
    already_calculated_q = Q(genecoveragecanonicaltranscript__isnull=False)
    #num_already_calculated = coverage_collection_qs.filter(already_calculated_q).distinct().count()
    #if num_already_calculated:
    #    logging.info("Skipping %d already calculated", num_already_calculated)

    for cc in coverage_collection_qs.exclude(already_calculated_q):
        transcript_ids, original_transcript_ids = canonical_transcripts
        qt = Q(transcript_id__in=transcript_ids)
        qrefseq = Q(original_transcript_id__in=original_transcript_ids)
        qs = cc.genecoverage_set.filter(qt | qrefseq)
        if qs.exists():
            #logging.info("Getting GeneCoverage records for %s", cc)
            canonical_transcripts_list = []

            for gc_dict in qs.values():
                # GeneCoverageCanonicalTranscript has all of GeneCoverage's fields
                del gc_dict['id']
                gc_dict["canonical_transcript_collection"] = canonical_collection
                canonical_coverage = GeneCoverageCanonicalTranscript(**gc_dict)
                canonical_transcripts_list.append(canonical_coverage)

            if canonical_transcripts_list:
                #logging.info("Bulk inserting %d GeneCoverageCanonicalTranscript records", len(canonical_transcripts_list))
                GeneCoverageCanonicalTranscript.objects.bulk_create(canonical_transcripts_list)
