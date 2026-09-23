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

What the file is imported for is the pair: its patient chain, the seqauto links and the specimen's
measures (@see upload.tso500.dragen_combined_variant_output_records). That hangs off a VCF and Sample
like every other arm file's, so DragenTSO500CombinedVariantOutputCreateVCFTask writes a VCF of no
records whose one sample is the RNA arm, and DragenTSO500CombinedVariantOutputInsertTask takes the rest
of the file once the header step has made the Sample. A chain that cannot be made is a message on the
import page rather than a failure.

The file names no genome build, so one is declared at upload (@see upload.upload_metadata) or comes
off the VCFSourceSettings row for '^DRAGEN TSO500 CombinedVariantOutput'.
"""
import logging

from library.genomics.vcf_writer import VCFWriter, build_header_lines
from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_LENGTH, GENE_LEVEL_CONTIG_NAME
from upload.models import SimpleVCFImportInfo, UploadStep
from upload.tasks.vcf.import_vcf_step_task import ImportVCFStepTask
from upload.tso500.dragen_combined_variant_output_parser import (
    MODULE_VERSION,
    RNA_SAMPLE_ID,
    get_analysis_details,
    read_combined_variant_output,
)
from upload.tso500.dragen_combined_variant_output_records import (
    CombinedVariantOutputIdentityError,
    link_samples_to_extractions,
    link_to_sequencing_run,
    measured_date,
    parse_pair_identifiers,
    resolve_pair,
    write_specimen_measures,
)
from variantgrid.celery import app

# The sample's FORMAT fields - read support rather than a genotype, bound by the
# '^DRAGEN TSO500 CombinedVariantOutput' VCFSourceSettings row. A header with no FORMAT has no sample
# columns, and the Sample would be named after the file rather than the RNA arm
ALT_READS_FORMAT = "ALT_READS"
REF_READS_FORMAT = "REF_READS"

# What '##source' says, so VCFSourceSettings can speak for this file type
SOURCE_PREFIX = "DRAGEN TSO500 CombinedVariantOutput"


def source_from_analysis_details(analysis_details: dict) -> str:
    """ 'DRAGEN TSO500 CombinedVariantOutput 2.1.1' - the software that wrote the file, which is
        what '##source' gives every other import """
    if module_version := analysis_details.get(MODULE_VERSION):
        return f"{SOURCE_PREFIX} {module_version}"
    return SOURCE_PREFIX


def _sample_name(analysis_details: dict, file_upload) -> str:
    """ The RNA arm, which is what links the pair to its sequencing run """
    return analysis_details.get(RNA_SAMPLE_ID) or file_upload.name


def _write_sample_vcf(filename: str, sample_name: str, source: str):
    """ A VCF of no records - only the Sample the rest of the file is written against """

    header_lines = build_header_lines(
        meta_lines=[f"##source={source}"] if source else [],
        formats=[
            f'##FORMAT=<ID={ALT_READS_FORMAT},Number=1,Type=Integer,'
            f'Description="Reads supporting the splice junction">',
            f'##FORMAT=<ID={REF_READS_FORMAT},Number=1,Type=Integer,'
            f'Description="Reads across the reference transcript at the junction">',
        ],
        contig_lines=[f"##contig=<ID={GENE_LEVEL_CONTIG_NAME},length={GENE_LEVEL_CONTIG_LENGTH}>"],
        samples=[sample_name],
    )
    with open(filename, "w") as f:
        VCFWriter(f, header_lines)


class DragenTSO500CombinedVariantOutputCreateVCFTask(ImportVCFStepTask):
    """ Write the VCF that makes the RNA arm's Sample, for the insert step to accession """

    def process_items(self, upload_step):
        sections = read_combined_variant_output(upload_step.input_filename)
        analysis_details = get_analysis_details(sections)
        file_upload = upload_step.upload_pipeline.file_upload
        _write_sample_vcf(upload_step.output_filename,
                          sample_name=_sample_name(analysis_details, file_upload),
                          source=source_from_analysis_details(analysis_details))
        return 0


class DragenTSO500CombinedVariantOutputInsertTask(ImportVCFStepTask):
    """ Runs once the VCF and its Sample exist - everything in the file that is not a variant: the
        pair's patient chain, the seqauto links and the specimen's measures """

    def process_items(self, upload_step: UploadStep):
        upload_pipeline = upload_step.upload_pipeline
        vcf = upload_pipeline.uploadedvcf.vcf
        file_upload = upload_pipeline.file_upload
        user = file_upload.user

        sections = read_combined_variant_output(file_upload.get_filename())
        analysis_details = get_analysis_details(sections)
        try:
            identifiers = parse_pair_identifiers(analysis_details)
            if identifiers is None:
                logging.info("%s names no pair to accession", file_upload)
                return 0
            resolved = resolve_pair(identifiers, user)
        except CombinedVariantOutputIdentityError as e:
            SimpleVCFImportInfo.add_message_count(1, str(e), upload_step)
            logging.warning("%s: %s", file_upload, e)
            return 0

        linked = link_samples_to_extractions(resolved, user)
        logging.info("%s: linked %d sample(s) to %s", file_upload, linked, resolved.specimen)

        # This VCF's one sample is the RNA arm's
        if identifiers.rna and (sample := vcf.sample_set.first()):
            if sequencing_run := link_to_sequencing_run(vcf, sample, identifiers.rna.sample_id):
                logging.info("%s: linked %s to %s", file_upload, sample, sequencing_run)
            else:
                message = f"No sequencing sample named '{identifiers.rna.sample_id}' - " \
                          f"VCF not linked to a sequencing run"
                SimpleVCFImportInfo.add_message_count(1, message, upload_step)

        measures = write_specimen_measures(sections, resolved, identifiers, user,
                                           method=source_from_analysis_details(analysis_details),
                                           date=measured_date(analysis_details))
        return len(measures)


DragenTSO500CombinedVariantOutputCreateVCFTask = app.register_task(
    DragenTSO500CombinedVariantOutputCreateVCFTask())
DragenTSO500CombinedVariantOutputInsertTask = app.register_task(
    DragenTSO500CombinedVariantOutputInsertTask())
