"""
The VCF writer the gene-level imports share (CNV, SpliceGirl, TSO 500 fusions): each resolves a caller's
output to GeneLevelRecords, and write_gene_level_vcf puts them on the gene-level contig with END = POS,
so they go through the normal insert pipeline (@see snpdb.gene_level_variants).
"""
from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any, Optional

import cyvcf2

from library.genomics.vcf_enums import VCFColumns
from library.genomics.vcf_writer import VCFInfoHeader, VCFWriter, build_header_lines
from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_LENGTH, GENE_LEVEL_CONTIG_NAME
from snpdb.vcf_utils import vcf_header_filter_ids

# Header lines the rewrite replaces with its own
_REPLACED_HEADER_PREFIXES = ("##fileformat", "##INFO", "##FORMAT", "##contig", "##ALT", "#CHROM")


@dataclass(kw_only=True)
class GeneLevelRecord:
    """ One rewritten record - the gene-level event it resolved to (anything with a variant_coordinate),
        the task's INFO, and the caller's columns kept as written """
    event: Any
    info: dict
    vcf_id: str = "."
    qual: str = "."
    vcf_filter: str = "."
    fmt: Optional[str] = None
    sample_calls: list[str] = field(default_factory=list)

    @staticmethod
    def from_caller(record: cyvcf2.Variant, event, info: dict) -> 'GeneLevelRecord':
        columns = str(record).rstrip("\n").split("\t")
        return GeneLevelRecord(event=event, info=info, vcf_id=columns[VCFColumns.ID],
                               qual=columns[VCFColumns.QUAL], vcf_filter=columns[VCFColumns.FILTER],
                               fmt=columns[VCFColumns.FORMAT], sample_calls=columns[VCFColumns.FORMAT + 1:])


def write_gene_level_vcf(filename: str, records: list[GeneLevelRecord], *, info: Iterable[VCFInfoHeader],
                         samples: Iterable[str], raw_header: str = "", meta_lines: Iterable[str] = (),
                         formats: Iterable[str] = (), encode_info=None):
    """ Written already-clean and sorted, which is what lets preprocess skip straight to the split.
        END = POS gives svlen 0 through vcf_get_ref_alt_svlen_and_modification, which needs one of
        SVLEN/END for any symbolic alt (and reads SVLEN=0 as absent).

        The caller's raw_header keeps its meta, contig and FORMAT lines, so the build is detected from the
        file as it would have been and VCFSourceSettings still binds its fields. A FILTER a record uses
        that nothing declares is declared here, as preprocess skips vcf_clean_and_filter, which would
        have (SpliceGirl writes LowUniqueAlignments undeclared) """

    caller_lines = raw_header.splitlines()
    meta_lines = [*meta_lines,
                  *(line for line in caller_lines if not line.startswith(_REPLACED_HEADER_PREFIXES))]
    used_filters = {f for record in records for f in record.vcf_filter.split(";")}
    for filter_id in sorted(used_filters - vcf_header_filter_ids(meta_lines) - {"PASS", "."}):
        meta_lines.append(f'##FILTER=<ID={filter_id},Description="Not declared by the caller">')

    header_lines = build_header_lines(
        meta_lines=meta_lines,
        info=[VCFInfoHeader(id="END", type="Integer", description="Stop position of the interval"), *info],
        formats=[*formats, *(line for line in caller_lines if line.startswith("##FORMAT"))],
        contig_lines=[*(line for line in caller_lines if line.startswith("##contig")),
                      f"##contig=<ID={GENE_LEVEL_CONTIG_NAME},length={GENE_LEVEL_CONTIG_LENGTH}>"],
        samples=samples,
    )

    with open(filename, "w") as f:
        writer = VCFWriter(f, header_lines, encode_info=encode_info)
        for record in sorted(records, key=lambda r: (r.event.variant_coordinate.position,
                                                     r.event.variant_coordinate.alt)):
            vc = record.event.variant_coordinate
            writer.write_record(vc.chrom, vc.position, vc.ref, vc.alt,
                                vcf_id=record.vcf_id, qual=record.qual, vcf_filter=record.vcf_filter,
                                info={"END": vc.position, **record.info}, fmt=record.fmt,
                                sample_calls=record.sample_calls)
