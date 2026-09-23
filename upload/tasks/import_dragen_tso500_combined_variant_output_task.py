"""
Loader for Illumina DRAGEN TSO 500's CombinedVariantOutput.tsv - one vendor's format, not a standard.
Anything here that reads a named section or column belongs to that format.

None of its sections is a variant source. Each is carried better elsewhere: small variants and copy
number on their own VCFs, fusions on AllFusions.csv, which keeps the caller, score, filters and
split/pair breakdown this file drops, and splice calls on the RNA arm's SpliceVariants.vcf
(@see upload.tasks.import_splicegirl_vcf_task).

'[Splice Variants]' is a filtered view of that VCF (SpliceGirl), row for row:
Breakpoint 1 = POS, Breakpoint 2 = INFO/END, Splice Supporting Reads = ALTDEDUP (FORMAT AD) and
Reference Reads Transcript = REFDEDUP (FORMAT DP), no off-by-one - confirmed on a real pair
(SACGF/variantgrid_sapath#457). What Illumina keeps is "passing splice variants that are contained
on genes EGFR, MET, and AR" (DRAGEN TSO 500 v2.5 Combined Variant Output,
https://help.tso500software.illumina.com/dragen-tso-500-guides/dragen-tso-500-v2.5/analysis-output/combined-variant-output):
a gene filter on top of FILTER=PASS, not a whitelist of junctions, so EGFRvII or vIVa qualify as
much as vIII, and a PASS call in any other panel gene is in the VCF but never here. Scientists differ
on whether they work from the filtered or unfiltered set, so the splice calls come from the VCF,
every record with its FILTER (#1903).

What the file is imported for is the pair: its patient chain, the seqauto links and the analysis
itself - one seqauto.models.DragenTSO500CombinedVariantOutput per (run, pair), holding its TMB, MSI and
GIS (@see upload.tso500.dragen_combined_variant_output_records).

Entry point is ImportDragenTSO500CombinedVariantOutputTask, a single-shot ImportTask (the file has no
variants and no coordinates, so there is no VCF pipeline) returning the number of rows written, as the
run's MetricsOutput does (@see upload.tasks.import_dragen_tso500_metrics_output_task).

The run is the upload's 'sequencing_run' metadata - the file names its run 'NA'. Without it the run
whose current sample sheet names the pair's sample IDs is used, and a pair no registered run names
cannot be keyed, so it fails the import saying so.

A chain that cannot be made is not a failure: the row is still the analysis', with its specimen claim
parked saying why, which the pair's page shows. reconcile_pending_extractions is fired afterwards so a
claim that only became resolvable with this file attaches straight away.
"""
import logging

from patients.tasks.extraction_matching_tasks import reconcile_pending_extractions
from seqauto.models import SequencingRun
from upload.models import UploadedDragenTSO500CombinedVariantOutput
from upload.tasks.import_task import ImportTask
from upload.tso500.dragen_combined_variant_output_parser import (
    DNA_SAMPLE_ID,
    PAIR_ID,
    RNA_SAMPLE_ID,
    get_analysis_details,
    read_combined_variant_output,
)
from upload.tso500.dragen_combined_variant_output_records import (
    CombinedVariantOutputIdentityError,
    link_samples_to_extractions,
    parse_pair_identifiers,
    resolve_pair,
    sequencing_run_for_sample_ids,
    write_combined_variant_output,
)
from upload.upload_metadata import SEQUENCING_RUN, get_metadata_sequencing_run_name
from variantgrid.celery import app


class ImportDragenTSO500CombinedVariantOutputTask(ImportTask):

    def process_items(self, file_upload):
        user = file_upload.user
        sections = read_combined_variant_output(file_upload.get_filename())
        analysis_details = get_analysis_details(sections)

        pair_id = analysis_details.get(PAIR_ID)
        if not pair_id:
            raise ValueError(f"{file_upload} names no 'Pair ID' - there is no pair to record the "
                             f"analysis against")

        sample_ids = [sample_id for sample_id in (analysis_details.get(DNA_SAMPLE_ID),
                                                  analysis_details.get(RNA_SAMPLE_ID)) if sample_id]
        if sequencing_run_name := get_metadata_sequencing_run_name(file_upload):
            sequencing_run = SequencingRun.objects.filter(name=sequencing_run_name).first()
        elif sequencing_run := sequencing_run_for_sample_ids(sample_ids):
            sequencing_run_name = sequencing_run.name
        else:
            raise ValueError(f"CombinedVariantOutput.tsv names its run 'NA' and no sequencing run "
                             f"names {', '.join(sample_ids)} - upload it with '{SEQUENCING_RUN}' "
                             f"metadata naming the run it came off")
        UploadedDragenTSO500CombinedVariantOutput.objects.get_or_create(file_upload=file_upload)

        identifiers = None
        resolved = None
        chain_error = None
        try:
            if identifiers := parse_pair_identifiers(analysis_details):
                resolved = resolve_pair(identifiers, user)
            else:
                chain_error = f"{file_upload} names no pair and specimen to accession"
                logging.info(chain_error)
        except CombinedVariantOutputIdentityError as e:
            chain_error = str(e)
            logging.warning("%s: %s", file_upload, e)

        if resolved:
            linked = link_samples_to_extractions(resolved, user)
            logging.info("%s: linked %d sample(s) to %s", file_upload, linked, resolved.specimen)

        cvo = write_combined_variant_output(sections, pair_id, sequencing_run, sequencing_run_name, user,
                                            resolved=resolved,
                                            specimen_reference=identifiers.specimen_reference if identifiers else None,
                                            parked_error=chain_error, file_upload=file_upload)
        logging.info("%s: recorded %s", file_upload, cvo)
        reconcile_pending_extractions.delay()
        return 1


ImportDragenTSO500CombinedVariantOutputTask = app.register_task(
    ImportDragenTSO500CombinedVariantOutputTask())  # @UndefinedVariable
