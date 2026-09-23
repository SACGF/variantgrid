"""
Loader for Illumina DRAGEN TSO 500's CombinedVariantOutput.tsv - one vendor's format, not a standard.
Anything here that reads a named section or column belongs to that format; the splice identity it
resolves to does not (@see genes.gene_splice).

Only the '[Splice Variants]' section is a variant source. It is what a scientist reports a splice call from,
and it names the gene and the two breakpoints rather than pretending to be a deletion. The other
sections are carried better elsewhere: small variants and copy number on their own VCFs, and fusions
on AllFusions.csv, which keeps the caller, score, filters and split/pair breakdown this file drops.

The section is a filtered view of the RNA arm's SpliceVariants.vcf (SpliceGirl), row for row:
Breakpoint 1 = POS, Breakpoint 2 = INFO/END, Splice Supporting Reads = ALTDEDUP (FORMAT AD) and
Reference Reads Transcript = REFDEDUP (FORMAT DP), no off-by-one - confirmed on a real pair
(SACGF/variantgrid_sapath#457). What Illumina keeps is "passing splice variants that are contained
on genes EGFR, MET, and AR" (DRAGEN TSO 500 v2.5 Combined Variant Output,
https://help.tso500software.illumina.com/dragen-tso-500-guides/dragen-tso-500-v2.5/analysis-output/combined-variant-output):
a gene filter on top of FILTER=PASS, not a whitelist of junctions, so EGFRvII or vIVa qualify as
much as vIII, and a PASS call in any other panel gene is in the VCF but never here. So far every
call SA Pathology has signed out is in one of the three (SACGF/variantgrid#1875), but nothing
stops SpliceGirl passing one elsewhere.

The rows become a VCF of gene-level variants which goes through the normal VCF import pipeline, so
the VCF/Sample/Cohort come from the header the way every other import's do, and the CohortGenotype
rows are written by the same SQL COPY path. Only the bcftools stages are skipped, since they all need
a reference base a gene-level locus does not have. @see snpdb.gene_level_variants for why these are
Variants at all, and upload.vcf.gene_level_vcf_preprocess for exactly what is skipped and why.

Two steps. DragenTSO500CombinedVariantOutputCreateVCFTask writes the VCF; a splice event has no
record of its own the way a fusion has a GeneFusion, so the alt carries the gene and the label and the
caller's row rides along in INFO. DragenTSO500CombinedVariantOutputInsertTask then takes the rest of
the file - the pair's patient chain, the seqauto links and the specimen's measures
(@see upload.tso500.dragen_combined_variant_output_records) - which needs the Sample, so it runs once
the header step has made it, alongside data insertion. Most pairs have no splice call, and their
file is still what accessions them. A chain that cannot be made is a message on the import page
rather than a failure: the splice calls are worth having whether or not the pair has been
accessioned yet.

The file names no genome build, and the create-VCF step needs one to place the breakpoints, so it is
declared at upload (@see upload.upload_metadata) or comes off the VCFSourceSettings row for
'^DRAGEN TSO500 CombinedVariantOutput'. That row also binds the sample's ALT_READS/REF_READS as alt
and ref depth, so the sample node's minimum-reads threshold and allele frequency work on splice calls;
the frequency it yields is the junction ratio.
"""
import logging

import simplejson

from genes.gene_fusions import parse_breakpoint
from genes.gene_splice import ResolvedSpliceEvent, SpliceEventResolver
from library.genomics.vcf_writer import (
    VCFInfoHeader,
    VCFWriter,
    build_header_lines,
    percent_encode_info_value,
)
from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_LENGTH, GENE_LEVEL_CONTIG_NAME
from snpdb.models import GenomeBuild
from upload.models import SimpleVCFImportInfo, UploadStep
from upload.tasks.vcf.import_vcf_step_task import ImportVCFStepTask
from upload.tso500.dragen_combined_variant_output_parser import (
    BREAKPOINT_1,
    BREAKPOINT_2,
    GENE,
    MODULE_VERSION,
    REFERENCE_READS_TRANSCRIPT,
    RNA_SAMPLE_ID,
    SPLICE_INFO,
    SPLICE_OBSERVATION_INFO,
    SPLICE_SUPPORTING_READS,
    get_analysis_details,
    get_splice_rows,
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
from upload.vcf.vcf_import import resolve_genome_build_from_source
from variantgrid.celery import app

# The sample's FORMAT fields - read support rather than a genotype
ALT_READS_FORMAT = "ALT_READS"
REF_READS_FORMAT = "REF_READS"
VCF_MISSING_VALUE = "."

# What '##source' says, so VCFSourceSettings can speak for this file type
SOURCE_PREFIX = "DRAGEN TSO500 CombinedVariantOutput"


def source_from_analysis_details(analysis_details: dict) -> str:
    """ 'DRAGEN TSO500 CombinedVariantOutput 2.1.1' - the software that wrote the file, which is
        what '##source' gives every other import """
    if module_version := analysis_details.get(MODULE_VERSION):
        return f"{SOURCE_PREFIX} {module_version}"
    return SOURCE_PREFIX


def _sample_name(analysis_details: dict, file_upload) -> str:
    """ The splice caller runs on the RNA arm, so its sample is the one the calls belong to """
    return analysis_details.get(RNA_SAMPLE_ID) or file_upload.name


def _read_count(observation: dict, column: str) -> str:
    try:
        return str(int(observation[column]))
    except (KeyError, TypeError, ValueError):
        return VCF_MISSING_VALUE


def _events_by_variant_coordinate(rows: list[dict], genome_build: GenomeBuild) -> dict:
    """ {resolved splice event: the caller's row}. Two rows resolving to one identity would be the
        same junction written twice, so the first stands """

    resolver = SpliceEventResolver(genome_build)
    events: dict[ResolvedSpliceEvent, dict] = {}
    for row in rows:
        breakpoints = [parse_breakpoint(row.get(c)) for c in (BREAKPOINT_1, BREAKPOINT_2)]
        if not all(breakpoints):
            logging.warning("Splice row %s has no usable breakpoints - skipped", row)
            continue
        (chrom, donor), (_acceptor_chrom, acceptor) = breakpoints
        resolved = resolver.resolve(row[GENE], chrom, donor, acceptor)
        if resolved is None:
            logging.warning("Splice row %s could not be resolved to a gene and contig - skipped", row)
            continue
        events.setdefault(resolved, row)
    return events


def _write_gene_level_vcf(filename: str, events: dict, sample_name: str, source: str):
    """ Written already-clean and sorted, which is what lets preprocess skip straight to the split.
        END = POS gives svlen 0 through vcf_get_ref_alt_svlen_and_modification, which needs one of
        SVLEN/END for any symbolic alt (and reads SVLEN=0 as absent). """

    header_lines = build_header_lines(
        meta_lines=[f"##source={source}"] if source else [],
        info=[
            VCFInfoHeader(id="END", type="Integer", description="Stop position of the interval"),
            VCFInfoHeader(id=SPLICE_INFO, type="String",
                          description="Splice event as gene and label, eg 'AR-V7'"),
            VCFInfoHeader(id=SPLICE_OBSERVATION_INFO, type="String",
                          description="JSON of the caller row this splice event was called from"),
        ],
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
        writer = VCFWriter(f, header_lines, encode_info=percent_encode_info_value)
        for resolved in sorted(events, key=lambda r: (r.gene.pk, r.alt)):
            observation = events[resolved]
            variant_coordinate = resolved.variant_coordinate
            info = {
                "END": variant_coordinate.position,
                SPLICE_INFO: resolved.canonical_str,
                SPLICE_OBSERVATION_INFO: simplejson.dumps(observation, ignore_nan=True),
            }
            sample_call = ":".join([_read_count(observation, SPLICE_SUPPORTING_READS),
                                    _read_count(observation, REFERENCE_READS_TRANSCRIPT)])
            writer.write_record(variant_coordinate.chrom, variant_coordinate.position,
                                variant_coordinate.ref, variant_coordinate.alt,
                                info=info, fmt=f"{ALT_READS_FORMAT}:{REF_READS_FORMAT}",
                                sample_calls=[sample_call])


class DragenTSO500CombinedVariantOutputCreateVCFTask(ImportVCFStepTask):
    """ Write the splice variants as a VCF so they go through the normal insert pipeline """

    def process_items(self, upload_step):
        sections = read_combined_variant_output(upload_step.input_filename)
        analysis_details = get_analysis_details(sections)
        source = source_from_analysis_details(analysis_details)
        file_upload = upload_step.upload_pipeline.file_upload
        # The VCF this step writes is what the build would normally be resolved from, so the source
        # line has to answer it here - @see upload.vcf.vcf_import.resolve_genome_build
        genome_build = resolve_genome_build_from_source(source, file_upload)
        if genome_build is None:
            raise ValueError(f"{upload_step.input_filename} declares no genome build - send one as "
                             f"upload metadata, as the breakpoints are positions in one")

        rows = get_splice_rows(sections)
        events = _events_by_variant_coordinate(rows, genome_build)
        _write_gene_level_vcf(upload_step.output_filename, events,
                              sample_name=_sample_name(analysis_details, file_upload), source=source)
        return len(rows)


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

        # The splice caller runs on the RNA arm, so that is the sample this VCF's calls came off
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
