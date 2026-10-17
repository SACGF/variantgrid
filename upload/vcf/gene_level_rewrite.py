"""
The gene-level layer of the VCF rewrite that the gene-level imports share (CNV, SpliceGirl, TSO 500
fusions): each turns a caller's output into records on the gene-level contig with END = POS, so they
go through the normal insert pipeline (@see snpdb.gene_level_variants, upload.vcf.gene_level_vcf_preprocess).

Each task keeps its own parsing and gene resolution and the two INFO fields it writes; this owns the rest:
  GeneLevelRecord          - one record: its resolved event plus the caller's columns kept as written
  gene_level_header_lines  - the header, from the caller's own where there is one
  write_gene_level_vcf     - the records, sorted, with END and the caller's INFO
  report_skip_counts       - the records each task skipped, as messages on the import page
"""
from collections import Counter
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field
from typing import Any, Optional

import cyvcf2

from library.genomics.vcf_enums import VCFColumns
from library.genomics.vcf_writer import VCFInfoHeader, VCFWriter, build_header_lines
from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_LENGTH, GENE_LEVEL_CONTIG_NAME
from snpdb.vcf_utils import vcf_header_filter_ids
from upload.models.models import SimpleVCFImportInfo, UploadStep

# Header lines the rewrite replaces with its own
_REPLACED_HEADER_PREFIXES = ("##fileformat", "##INFO", "##FORMAT", "##contig", "##ALT", "#CHROM")

GENE_LEVEL_CONTIG_LINE = f"##contig=<ID={GENE_LEVEL_CONTIG_NAME},length={GENE_LEVEL_CONTIG_LENGTH}>"
END_INFO_HEADER = VCFInfoHeader(id="END", type="Integer", description="Stop position of the interval")


@dataclass(kw_only=True)
class GeneLevelRecord:
    """ One rewritten record - the gene-level event it resolved to (anything with a variant_coordinate),
        plus the caller's columns kept as written """
    event: Any
    vcf_id: str = "."
    qual: str = "."
    vcf_filter: str = "."
    fmt: Optional[str] = None
    sample_calls: list[str] = field(default_factory=list)

    @staticmethod
    def caller_columns(record: cyvcf2.Variant) -> dict:
        """ The caller's ID, QUAL, FILTER and sample columns as written, as GeneLevelRecord kwargs """
        columns = str(record).rstrip("\n").split("\t")
        return {
            "vcf_id": columns[VCFColumns.ID],
            "qual": columns[VCFColumns.QUAL],
            "vcf_filter": columns[VCFColumns.FILTER],
            "fmt": columns[VCFColumns.FORMAT],
            "sample_calls": columns[VCFColumns.FORMAT + 1:],
        }

    @property
    def sort_key(self) -> tuple:
        return self.event.gene.pk, self.event.alt


def _header_lines_from(raw_header: str, prefix: str) -> list[str]:
    return [line for line in raw_header.splitlines() if line.startswith(prefix)]


def gene_level_header_lines(*, info: Iterable[VCFInfoHeader], samples: Iterable[str],
                            records: Iterable[GeneLevelRecord], raw_header: Optional[str] = None,
                            meta_lines: Optional[Iterable[str]] = None,
                            formats: Optional[Iterable[str]] = None) -> list[str]:
    """ END plus the task's own INFO, on the gene-level contig.

        Given the caller's raw_header, its meta lines (source, FILTER...), contigs and FORMAT lines are
        kept, so the build is detected from the file as it would have been and VCFSourceSettings still
        binds its fields. A FILTER a record uses that nothing declares is declared here: the gene-level
        preprocess skips vcf_clean_and_filter, which would otherwise have carried it through (SpliceGirl
        writes LowUniqueAlignments undeclared). htslib adds a 'Dummy' declaration to a reader's raw_header
        as it reads such a record, so only a header read before the records is missing one """

    meta_lines = list(meta_lines or [])
    contig_lines = []
    if raw_header is not None:
        meta_lines = [line for line in raw_header.splitlines()
                      if not line.startswith(_REPLACED_HEADER_PREFIXES)] + meta_lines
        contig_lines = _header_lines_from(raw_header, "##contig")
        if formats is None:
            formats = _header_lines_from(raw_header, "##FORMAT")
    contig_lines.append(GENE_LEVEL_CONTIG_LINE)

    declared_filters = vcf_header_filter_ids(meta_lines)
    used_filters = {f for record in records for f in record.vcf_filter.split(";")}
    for filter_id in sorted(used_filters - declared_filters - {"PASS", "."}):
        meta_lines.append(f'##FILTER=<ID={filter_id},Description="Not declared by the caller">')

    return build_header_lines(meta_lines=meta_lines, info=[END_INFO_HEADER, *info], formats=formats,
                              contig_lines=contig_lines, samples=samples)


def write_gene_level_vcf(filename: str, header_lines: list[str], records: Iterable[GeneLevelRecord],
                         info_func: Callable[[GeneLevelRecord], dict], encode_info=None):
    """ Written already-clean and sorted, which is what lets preprocess skip straight to the split.
        END = POS gives svlen 0 through vcf_get_ref_alt_svlen_and_modification, which needs one of
        SVLEN/END for any symbolic alt (and reads SVLEN=0 as absent).

        info_func gives the task's own INFO for a record, written after END """

    with open(filename, "w") as f:
        writer = VCFWriter(f, header_lines, encode_info=encode_info)
        for record in sorted(records, key=lambda r: r.sort_key):
            variant_coordinate = record.event.variant_coordinate
            info = {"END": variant_coordinate.position, **info_func(record)}
            writer.write_record(variant_coordinate.chrom, variant_coordinate.position,
                                variant_coordinate.ref, variant_coordinate.alt,
                                vcf_id=record.vcf_id, qual=record.qual, vcf_filter=record.vcf_filter,
                                info=info, fmt=record.fmt, sample_calls=record.sample_calls)


def report_skip_counts(upload_step: UploadStep, skip_counts: Iterable[tuple[int, str]] = (), *,
                       ambiguous: Optional[Counter[str]] = None, ambiguous_prefix: str = ""):
    """ Each (count, message) with records skipped, then each AmbiguousGeneSymbolAlias message (which
        names the genes it could be) under ambiguous_prefix """
    for count, message in skip_counts:
        if count:
            SimpleVCFImportInfo.add_message_count(count, message, upload_step)
    for ambiguous_message, count in (ambiguous or {}).items():
        SimpleVCFImportInfo.add_message_count(count, ambiguous_prefix + ambiguous_message, upload_step)
