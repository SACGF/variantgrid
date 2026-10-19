"""
Loader for a CNV caller's VCF whose records are whole-gene calls - DRAGEN's TSO 500 CNV output, where
every record carries the gene it is about in INFO/SEGID.

What makes this a gene-level file is that segment field (settings.VCF_GENE_LEVEL_SEGMENT_FIELDS), not
the alt: the caller's gene claim is the signal. The segment coordinates it writes are the panel's
target window rather than the event - PIK3CD's segment starts 59 kb inside the gene, at the first
targeted exon - so they are not what the Variant is, @see genes.gene_copy_number for why identity is
the gene plus the direction. A caller that names a gene on partial, exon-level calls (DragenExonCNV's
GENE=BRCA1) keeps importing as coordinate structural variants; only the segment field claims a file.

The records are rewritten onto the gene-level contig and go through the normal VCF import pipeline,
so the VCF/Sample/Cohort come from the header the way every other import's do, and the CohortGenotype
rows are written by the same SQL COPY path. Only the bcftools stages are skipped, since they all need
a reference base a gene-level locus does not have (@see upload.vcf.gene_level_vcf_preprocess).

The rewrite keeps the caller's own header - its contigs (so the genome build is detected from the
file exactly as it would have been), its FILTER and FORMAT lines, and each record's ID, QUAL, FILTER
and sample columns as written. That is what leaves SM in FORMAT for
upload.vcf.vcf_import.get_copy_number_field to bind as the copy number column. What the record loses
is the caller's INFO, which is all segment geometry.

Two steps:
  GeneLevelCNVCreateVCFTask - rewrite the caller's VCF as gene-level records
  GeneLevelCNVInsertTask    - once the variants exist, create the GeneCopyNumberEvents against them
"""
import logging
from collections import Counter
from dataclasses import dataclass, field
from typing import Optional

import cyvcf2
from django.conf import settings

from genes.gene_copy_number import (
    ResolvedGeneCopyNumberEvent,
    create_gene_copy_number_events_for_variants,
)
from genes.gene_level_resolver import GeneLevelNameResolver
from genes.models import GeneCopyNumberEventKind
from genes.models.models_gene import AmbiguousGeneSymbolAlias
from library.genomics.vcf_enums import VCFSymbolicAllele
from library.genomics.vcf_writer import VCFInfoHeader
from upload.models import SimpleVCFImportInfo, UploadStep
from upload.tasks.vcf.import_vcf_step_task import ImportVCFStepTask
from upload.vcf.gene_level_rewrite import GeneLevelRecord, write_gene_level_vcf
from variantgrid.celery import app

# The canonical string of the event a record became, so the file says what it is - the peer of the
# fusion loader's FUSION. The segment field is written through with the gene name as the caller wrote
# it, which is what binds VCF.gene_level_segment_field when the header is read back
GENE_COPY_NUMBER_INFO = "GENE_CN"

# Which way a called record went. Only these two alts are a whole-gene copy number call
COPY_NUMBER_KIND_BY_ALT = {
    VCFSymbolicAllele.DUP: GeneCopyNumberEventKind.GAIN,
    VCFSymbolicAllele.DEL: GeneCopyNumberEventKind.LOSS,
}

VCF_MISSING_VALUE = "."

NO_SEGMENT_MESSAGE = "Gene-level CNV records skipped - no gene named in the segment field"
DUPLICATE_MESSAGE = "Gene-level CNV records skipped - the gene and direction were already called"
AMBIGUOUS_MESSAGE_PREFIX = "Gene-level CNV records skipped - "


def can_process_file(filename: str) -> bool:
    """ A VCF whose header declares one of the segment fields is a file of whole-gene calls """
    return bool(get_segment_field(filename))


def get_segment_field(filename: str) -> Optional[str]:
    """ Which of settings.VCF_GENE_LEVEL_SEGMENT_FIELDS this file names its genes in """
    try:
        reader = cyvcf2.VCF(filename)
    except Exception:  # Every VCF factory is asked, so a file we can't read is someone else's
        return None
    for field_id in settings.VCF_GENE_LEVEL_SEGMENT_FIELDS:
        if f"##INFO=<ID={field_id}," in reader.raw_header:
            return field_id
    return None


@dataclass
class GeneLevelCNVRewrite:
    """ What the rewrite made of the caller's file """
    records: list[GeneLevelRecord] = field(default_factory=list)
    records_read: int = 0
    no_segment: int = 0
    duplicates: int = 0
    ambiguous: Counter[str] = field(default_factory=Counter)  # AmbiguousGeneSymbolAlias message -> records


def _segment_gene(record: cyvcf2.Variant, segment_field: str) -> Optional[str]:
    """ The gene the caller's segment named, or None where it named none - VCF's own missing value
        included, since '.' is a spelling of blank rather than a gene """
    value = (record.INFO.get(segment_field) or "").strip()
    if not value or value == VCF_MISSING_VALUE:
        return None
    return value


def read_gene_level_cnv_records(reader: cyvcf2.VCF, segment_field: str,
                                resolver: GeneLevelNameResolver = None) -> GeneLevelCNVRewrite:
    """ The called <DUP>/<DEL> records that name a gene, as the gene-level events they will be
        stored as. A no-call is not an event, and a called record naming no gene is counted rather
        than guessed at, so a file this rule does not fit says so instead of half-importing.

        A name HGNC doesn't know is no reason to drop a call - it gets a custom GeneLevelId, as an
        unknown fusion partner does. An ambiguous alias is dropped and counted under its message, which
        names the genes it could be. """

    if resolver is None:
        resolver = GeneLevelNameResolver()

    rewrite = GeneLevelCNVRewrite()
    seen = set()
    for record in reader:
        rewrite.records_read += 1
        kind = COPY_NUMBER_KIND_BY_ALT.get(record.ALT[0]) if record.ALT else None
        if kind is None:
            continue  # A no-call, or an alt that isn't a whole-gene copy number call

        written_gene = _segment_gene(record, segment_field)
        try:
            resolved_gene = resolver.resolve_gene(written_gene) if written_gene else None
        except AmbiguousGeneSymbolAlias as ambiguous:
            rewrite.ambiguous[str(ambiguous)] += 1
            continue
        if resolved_gene is None:
            rewrite.no_segment += 1
            continue

        event = ResolvedGeneCopyNumberEvent(gene=resolved_gene.gene_level_id, kind=kind)
        if event in seen:
            # Two of the caller's segments resolved to one gene - MYCL1 and MYCL are the same gene,
            # and one Variant can carry one call
            rewrite.duplicates += 1
            continue
        seen.add(event)

        info = {segment_field: written_gene, GENE_COPY_NUMBER_INFO: event.canonical_str}
        rewrite.records.append(GeneLevelRecord.from_caller(record, event, info))
    return rewrite


def write_gene_level_cnv_vcf(filename: str, rewrite: GeneLevelCNVRewrite, reader: cyvcf2.VCF,
                             segment_field: str):
    write_gene_level_vcf(filename, rewrite.records, info=[
        VCFInfoHeader(id=segment_field, type="String", description="Gene the caller's segment named, as written"),
        VCFInfoHeader(id=GENE_COPY_NUMBER_INFO, type="String", description="Whole-gene copy number event"),
    ], samples=reader.samples, raw_header=reader.raw_header)


class GeneLevelCNVCreateVCFTask(ImportVCFStepTask):
    """ Rewrite the caller's whole-gene calls onto the gene-level contig, so they go through the
        normal insert pipeline """

    def process_items(self, upload_step: UploadStep):
        filename = upload_step.input_filename
        segment_field = get_segment_field(filename)
        if segment_field is None:
            raise ValueError(f"'{filename}' declares none of {settings.VCF_GENE_LEVEL_SEGMENT_FIELDS}")

        reader = cyvcf2.VCF(filename)
        rewrite = read_gene_level_cnv_records(reader, segment_field)
        write_gene_level_cnv_vcf(upload_step.output_filename, rewrite, reader, segment_field)

        for count, message in [(rewrite.no_segment, NO_SEGMENT_MESSAGE),
                               (rewrite.duplicates, DUPLICATE_MESSAGE)]:
            if count:
                SimpleVCFImportInfo.add_message_count(count, message, upload_step)
        for ambiguous_message, count in rewrite.ambiguous.items():
            SimpleVCFImportInfo.add_message_count(count, AMBIGUOUS_MESSAGE_PREFIX + ambiguous_message, upload_step)
        return rewrite.records_read


class GeneLevelCNVInsertTask(ImportVCFStepTask):
    """ Runs after data insertion, so the Variants and their CohortGenotypes exist - create the
        GeneCopyNumberEvents against them """

    def process_items(self, upload_step: UploadStep):
        vcf = upload_step.upload_pipeline.uploadedvcf.vcf
        created = create_gene_copy_number_events_for_variants(vcf.get_variant_qs())
        logging.info("Created %d gene copy number events for %s", created, vcf)
        return created


GeneLevelCNVCreateVCFTask = app.register_task(GeneLevelCNVCreateVCFTask())
GeneLevelCNVInsertTask = app.register_task(GeneLevelCNVInsertTask())
