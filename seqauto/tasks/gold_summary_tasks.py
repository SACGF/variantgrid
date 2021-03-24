"""

Based on code from TAU project - archive.make_gold_ref_standard

We want all coverage data, not just that which was properly matched converting between RefSeq/Ensembl
Do everything using originally provided gene_symbol / refseq_gene_id
Sometimes you can match the same RefSeq symbol to Ensembl - eg:

eg: MIR4697HG/NR_024344 and MIR4697/NR_039846 both => ENSG00000264919, ENST00000582977
RefSeq 'ANKRD18A' and 'FAM95C' both resolve to ENSG00000273170

"""
import celery
from collections import defaultdict
import logging
from scipy import stats

from eventlog.models import create_event
from genes.models import GeneCoverageCollection, GeneCoverageCanonicalTranscript
from library.enums.log_level import LogLevel
from library.log_utils import get_traceback, log_traceback
from library.utils import get_single_element
import numpy as np
from seqauto.models import SequencingRun, GoldReference, GoldCoverageSummary, GoldGeneCoverageCollection, EnrichmentKit
from snpdb.models import DataState
from snpdb.models.models_enums import ImportStatus


@celery.task
def calculate_gold_summary(enrichment_kit_id):
    """ Set SequencingRun.gold_standard on runs before this task is called """

    enrichment_kit = EnrichmentKit.objects.get(pk=enrichment_kit_id)
    logging.info("calculate_gold_summary for %s", enrichment_kit)
    GoldReference.objects.filter(enrichment_kit=enrichment_kit).delete()
    gold_reference = GoldReference.objects.create(enrichment_kit=enrichment_kit,
                                                  import_status=ImportStatus.IMPORTING)

    try:
        seq_runs_qs = SequencingRun.objects.filter(samplesheet__sequencingsample__enrichment_kit=enrichment_kit,
                                                   gold_standard=True).distinct()
        # Check that we have sufficient data from all of them...
        problems = []
        for sequencing_run in seq_runs_qs:
            try:
                sample_sheet = sequencing_run.sequencingruncurrentsamplesheet.sample_sheet
                for ss in sample_sheet.sequencingsample_set.filter(is_control=False):
                    qc = ss.get_single_qc()
                    if qc is None:
                        problems.append(f"{sequencing_run} - {ss} has no QC")
                    else:
                        try:
                            qc_gene_coverage = qc.qcgenecoverage
                            if qc_gene_coverage.data_state != DataState.COMPLETE:
                                msg = f"{sequencing_run} - {ss} QCGeneCoverage status = {qc_gene_coverage.get_data_state_display()}"
                                problems.append(msg)
                                continue

                            # Track what's used
                            GoldGeneCoverageCollection.objects.create(gold_reference=gold_reference,
                                                                      gene_coverage_collection=qc_gene_coverage.gene_coverage_collection)
                        except GeneCoverageCollection.DoesNotExist:
                            problems.append(f"{sequencing_run} - {ss} has no GeneCoverageCollection")

            except:
                problems.append(get_traceback())

        if problems:
            for p in problems:
                logging.error(p)
            error_exception = '\n'.join(problems)
            gold_reference.error_exception = error_exception
            gold_reference.import_status = ImportStatus.ERROR
            create_event(None, "calculate_gold_summary", error_exception, severity=LogLevel.ERROR)
        else:
            seq_run_ids = seq_runs_qs.values_list("pk", flat=True)
            gold_coverage_summaries = get_gold_coverage_summaries(gold_reference, seq_run_ids)
            GoldCoverageSummary.objects.bulk_create(gold_coverage_summaries)
            gold_reference.import_status = ImportStatus.SUCCESS
    except:
        error_exception = get_traceback()
        details = "enrichment_kit_id=%d\n" % enrichment_kit_id
        details += error_exception
        create_event(None, "calculate_gold_summary", details, severity=LogLevel.ERROR)
        gold_reference.error_exception = error_exception
        gold_reference.import_status = ImportStatus.ERROR

    gold_reference.save()


def get_gold_coverage_summaries(gold_reference, seq_run_ids):
    seq_run_path = "gene_coverage_collection__qcgenecoverage__qc__bam_file__unaligned_reads"
    seq_run_path += "__sequencing_sample__sample_sheet__sequencingruncurrentsamplesheet__sequencing_run"
    kwargs = {seq_run_path + "__in": seq_run_ids}
    qs = GeneCoverageCanonicalTranscript.objects.filter(**kwargs)

    # If there was a transcript ID match, there will be exactly 1 gene
    # Otherwise, there may be multiple genes with no transcript_id
    gene_symbol_transcript_arrays = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    for data in qs.values():
        original_gene_symbol = data["original_gene_symbol"]
        original_transcript_id = data["original_transcript_id"]
        arrays = gene_symbol_transcript_arrays[original_gene_symbol][original_transcript_id]

        for k, v in data.items():
            arrays[k].append(v)

    return calculate_gold_summary_stats(gold_reference, gene_symbol_transcript_arrays)


def calculate_gold_summary_stats(gold_reference, gene_symbol_transcript_gene_arrays):
    gold_coverage_summaries = []
    for original_gene_symbol, transcript_arrays in gene_symbol_transcript_gene_arrays.items():
        for original_transcript_id, arrays_dict in transcript_arrays.items():
            min_array = np.array(arrays_dict["min"])
            mean_array = np.array(arrays_dict["mean"])
            percent_10x_array = np.array(arrays_dict["percent_10x"])
            percent_20x_array = np.array(arrays_dict["percent_20x"])

            #5th percentile
            depth_10x_5th_percentile = np.percentile(percent_10x_array, 5)
            depth_20x_5th_percentile = np.percentile(percent_20x_array, 5)
            depth_mean_5th_percentile = np.percentile(mean_array, 5)
            depth_mean_95th_percentile = np.percentile(mean_array, 95)

            #get mean/std of mins
            min_mean = np.mean(min_array)
            #min_depth_std = np.std(min_depth_array)

            #get mean/std of means
            mean = np.mean(mean_array)
            #mean_depth_std = np.std(mean_array)

            standard_error = stats.sem(mean_array)

            try:
                # These should have all matched consistently...
                gene_symbol = get_single_element(set(arrays_dict["gene_symbol"]))
                transcript_id = get_single_element(set(arrays_dict["transcript_id"]))
            except:
                log_traceback()
                logging.error(arrays_dict)
                raise

            gcs = GoldCoverageSummary(gold_reference=gold_reference,
                                      gene_symbol=gene_symbol,
                                      transcript_id=transcript_id,
                                      original_gene_symbol=original_gene_symbol,
                                      original_transcript_id=original_transcript_id,
                                      mean=mean,
                                      standard_error=standard_error,
                                      min_mean=min_mean,
                                      depth_20x_5th_percentile=depth_20x_5th_percentile,
                                      depth_10x_5th_percentile=depth_10x_5th_percentile,
                                      depth_mean_5th_percentile=depth_mean_5th_percentile,
                                      depth_mean_95th_percentile=depth_mean_95th_percentile)
            gold_coverage_summaries.append(gcs)

    return gold_coverage_summaries
