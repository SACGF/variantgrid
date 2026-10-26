"""
Loader for SpliceGirl's SpliceVariants.vcf - the DRAGEN TSO 500 RNA arm's splice caller - as gene-level
splice variants (@see genes.gene_splice, and snpdb.gene_level_variants for why a splice call is a
Variant at all).

SpliceGirl writes each junction as a <DEL> from POS (Breakpoint 1) to INFO/END (Breakpoint 2). Imported
as it stands that is a genomic deletion: annotated as one, with no splice label, and on a Locus whose
REF is the caller's base at POS+2 rather than the reference base (upload/test_data/tso500/README.md).
What the file is recognised by is its '##source=SpliceGirl' header line, the same thing the
'^SpliceGirl' VCFSourceSettings row matches - the pipeline sends it as a plain VCF.

Every record is kept with its FILTER, not just the PASS calls: scientists differ on whether they work
from the filtered or unfiltered set, so the analysis filters decide. The CombinedVariantOutput's
'[Splice Variants]' is the PASS records on EGFR, MET and AR, and is no longer a variant source
(@see seqauto.tasks.import_dragen_tso500_combined_variant_output_task).

The record names no gene. It comes from the SpliceEvent row at the junction's coordinates, else the
one gene a transcript of the build puts there (@see genes.gene_splice.SpliceEventResolver.resolve_junction).

The rewrite keeps the caller's own header - its contigs (so the build is detected from the file as
it would have been), its source, QUAL, FILTER and FORMAT lines - and each record's ID, QUAL, FILTER and
sample columns as written. So the '^SpliceGirl' VCFSourceSettings row still binds AD/DP as alt/ref
depth, and VAF is the junction ratio ALTDEDUP / (ALTDEDUP + REFDEDUP). The caller's INFO rides along in
the splice observation, written with the CombinedVariantOutput's column names so the analysis grid's
'Splice calls' column reads either.

Only the bcftools stages are skipped, since they need a reference base a gene-level locus does not
have (@see upload.vcf.gene_level_vcf_preprocess).
"""
from dataclasses import dataclass, field
from typing import Optional

import cyvcf2

from genes.gene_splice import ResolvedSpliceEvent, SpliceEventResolver
from library.genomics.vcf_enums import VCFSymbolicAllele
from library.genomics.vcf_writer import VCFInfoHeader, percent_encode_info_value
from library.utils.json_utils import json_dumps_nan_as_null
from snpdb.models import GenomeBuild
from upload.models import SimpleVCFImportInfo, UploadStep
from upload.tasks.vcf.import_vcf_step_task import ImportVCFStepTask
from upload.tso500.dragen_combined_variant_output_parser import (
    BREAKPOINT_1,
    BREAKPOINT_2,
    GENE,
    REFERENCE_READS_TRANSCRIPT,
    SPLICE_INFO,
    SPLICE_OBSERVATION_INFO,
    SPLICE_SUPPORTING_READS,
)
from upload.vcf.gene_level_rewrite import GeneLevelRecord, write_gene_level_vcf
from upload.vcf.vcf_import import resolve_genome_build
from variantgrid.celery import app

# What '##source' starts with - the same text the '^SpliceGirl' VCFSourceSettings row matches
SPLICEGIRL_SOURCE_PREFIX = "SpliceGirl"
SPLICE_SUPPORTING_READS_INFO = "ALTDEDUP"
REFERENCE_READS_INFO = "REFDEDUP"

NOT_A_JUNCTION_MESSAGE = "SpliceGirl records skipped - not a <DEL> with an END"
UNRESOLVED_MESSAGE = "SpliceGirl records skipped - no single gene at the junction"
DUPLICATE_MESSAGE = "SpliceGirl records skipped - the junction was already called"


def can_process_file(filename: str) -> bool:
    try:
        reader = cyvcf2.VCF(filename)
    except Exception:  # Every VCF factory is asked, so a file we can't read is someone else's
        return False
    return any(line.startswith(f"##source={SPLICEGIRL_SOURCE_PREFIX}")
               for line in reader.raw_header.splitlines())


@dataclass
class SpliceGirlRewrite:
    """ What the rewrite made of the caller's file """
    records: list[GeneLevelRecord] = field(default_factory=list)
    records_read: int = 0
    not_a_junction: int = 0
    unresolved: int = 0
    duplicates: int = 0


def _junction_end(record: cyvcf2.Variant) -> Optional[int]:
    if record.ALT and record.ALT[0] == VCFSymbolicAllele.DEL:
        if end := record.INFO.get("END"):
            return int(end)
    return None


def _observation(record: cyvcf2.Variant, event: ResolvedSpliceEvent, end: int) -> dict:
    """ The caller's record under the CombinedVariantOutput's column names, which is what
        format_splice_observation renders, then the caller's INFO as written """
    observation = {
        GENE: event.gene.symbol_str,
        BREAKPOINT_1: f"{record.CHROM}:{record.POS}",
        BREAKPOINT_2: f"{record.CHROM}:{end}",
        SPLICE_SUPPORTING_READS: record.INFO.get(SPLICE_SUPPORTING_READS_INFO),
        REFERENCE_READS_TRANSCRIPT: record.INFO.get(REFERENCE_READS_INFO),
    }
    observation.update(dict(record.INFO))
    return observation


def read_splicegirl_records(reader: cyvcf2.VCF, genome_build: GenomeBuild,
                            resolver: SpliceEventResolver = None) -> SpliceGirlRewrite:
    """ Every junction record as the gene-level splice event it will be stored as. A record that
        resolves to no gene is counted rather than guessed at, so the import page says so """

    if resolver is None:
        resolver = SpliceEventResolver(genome_build)

    rewrite = SpliceGirlRewrite()
    seen = set()
    for record in reader:
        rewrite.records_read += 1
        if (end := _junction_end(record)) is None:
            rewrite.not_a_junction += 1
            continue

        event = resolver.resolve_junction(record.CHROM, record.POS, end)
        if event is None:
            rewrite.unresolved += 1
            continue
        if event in seen:
            rewrite.duplicates += 1
            continue
        seen.add(event)

        info = {SPLICE_INFO: event.canonical_str,
                SPLICE_OBSERVATION_INFO: json_dumps_nan_as_null(_observation(record, event, end))}
        rewrite.records.append(GeneLevelRecord.from_caller(record, event, info))
    return rewrite


def write_splicegirl_vcf(filename: str, rewrite: SpliceGirlRewrite, reader: cyvcf2.VCF):
    write_gene_level_vcf(filename, rewrite.records, info=[
        VCFInfoHeader(id=SPLICE_INFO, type="String", description="Splice event as gene and label, eg 'AR-V7'"),
        VCFInfoHeader(id=SPLICE_OBSERVATION_INFO, type="String",
                      description="JSON of the caller record this splice event was called from"),
    ], samples=reader.samples, raw_header=reader.raw_header, encode_info=percent_encode_info_value)


class SpliceGirlCreateVCFTask(ImportVCFStepTask):
    """ Rewrite SpliceGirl's junctions onto the gene-level contig, so they go through the normal
        insert pipeline """

    def process_items(self, upload_step: UploadStep):
        filename = upload_step.input_filename
        reader = cyvcf2.VCF(filename)
        genome_build = resolve_genome_build(reader, upload_step.upload_pipeline.file_upload)
        if genome_build is None:
            raise ValueError(f"{filename} declares no genome build - send one as upload metadata, "
                             f"as the breakpoints are positions in one")

        rewrite = read_splicegirl_records(reader, genome_build)
        write_splicegirl_vcf(upload_step.output_filename, rewrite, cyvcf2.VCF(filename))

        for count, message in [(rewrite.not_a_junction, NOT_A_JUNCTION_MESSAGE),
                               (rewrite.unresolved, UNRESOLVED_MESSAGE),
                               (rewrite.duplicates, DUPLICATE_MESSAGE)]:
            if count:
                SimpleVCFImportInfo.add_message_count(count, message, upload_step)
        return rewrite.records_read


SpliceGirlCreateVCFTask = app.register_task(SpliceGirlCreateVCFTask())
