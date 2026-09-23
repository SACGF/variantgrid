"""
Reader for Illumina DRAGEN TSO 500's CombinedVariantOutput.tsv - the file only, no database work.

What becomes of it is upload.tasks.import_dragen_tso500_combined_variant_output_task. The
'[Splice Variants]' rows are the only variants read from here - small variants and copy number stay on
their VCFs and fusions on AllFusions.csv, which carries the caller, score and filters this file drops -
and '[Analysis Details]', '[TMB]', '[MSI]' and '[GIS]' name the pair and its measures
(@see upload.tso500.dragen_combined_variant_output_records).

The file is the pair-level summary written beside the two arm directories: a banner, then one
'[Section]' per call type. A section is either key/value ('[Analysis Details]', '[TMB]') or a table
whose first line is its header, and it ends at a blank line. Every line is right-padded with tabs to
the widest section, and a lone 'NA' stands for a section with nothing in it.

Section names move between module versions - 2.1.1's '[Exon-Level CNVs]' is documented as
'Large Rearrangements' in 2.6, which also adds 'Gene-level Loss of Heterozygosity' - so a section
this has no name for is read and kept rather than being an error.
"""
import re
from dataclasses import dataclass, field
from typing import Optional

# The banner line, which is what recognises the file - a tsv of gene symbols would otherwise be
# claimed as a gene list. 2.6 puts its version in it, eg 'DRAGEN TruSight Oncology 500 v2.6.2 Analysis...'
FIRST_LINE_PATTERN = re.compile(r"DRAGEN TruSight Oncology 500 (?:v[\d.]+ )?Analysis Software - Combined Variant Output")

ANALYSIS_DETAILS = "Analysis Details"
SPLICE_VARIANTS = "Splice Variants"
TMB = "TMB"
MSI = "MSI"
GIS = "GIS"

# [Analysis Details] keys
PAIR_ID = "Pair ID"
DNA_SAMPLE_ID = "DNA Sample ID"
RNA_SAMPLE_ID = "RNA Sample ID"
OUTPUT_DATE = "Output Date"
OUTPUT_TIME = "Output Time"
MODULE_VERSION = "Module Version"
PIPELINE_VERSION = "Pipeline Version"

# The key/value sections holding the pair's measures - what each is written as is
# upload.tso500.dragen_combined_variant_output_records
TOTAL_TMB = "Total TMB"
CODING_REGION_SIZE = "Coding Region Size in Megabases"
PASSING_ELIGIBLE_VARIANTS = "Number of Passing Eligible Variants"
USABLE_MSI_SITES = "Usable MSI Sites"
TOTAL_MSI_SITES_UNSTABLE = "Total MSI Sites Unstable"
PERCENT_UNSTABLE_MSI_SITES = "Percent Unstable MSI Sites"
GENOMIC_INSTABILITY_SCORE = "Genomic Instability Score"
TUMOR_FRACTION = "Tumor Fraction"
PLOIDY = "Ploidy"

# [Splice Variants] columns
GENE = "Gene"
AFFECTED_EXON = "Affected Exon"
BREAKPOINT_1 = "Breakpoint 1"
BREAKPOINT_2 = "Breakpoint 2"
SPLICE_SUPPORTING_READS = "Splice Supporting Reads"
REFERENCE_READS_TRANSCRIPT = "Reference Reads Transcript"

SPLICE_COLUMNS = (GENE, BREAKPOINT_1, BREAKPOINT_2)

# INFO fields a splice row is carried in once it becomes a VCF - they land in CohortGenotype.info via
# the standard bulk importer, which stores every INFO field the header declares
SPLICE_INFO = "SPLICE"
SPLICE_OBSERVATION_INFO = "SPLICE_OBS"

MISSING_VALUE = "NA"


@dataclass
class CombinedVariantOutputSection:
    """ One '[Section]' of the file. A key/value section reads `values`, a table reads
        `header`/`rows`; `lines` is what both are derived from, so a section we have no name for
        still arrives whole """
    name: str
    lines: list[list[str]] = field(default_factory=list)

    @property
    def header(self) -> list[str]:
        return self.lines[0] if self.lines else []

    @property
    def rows(self) -> list[dict]:
        """ The table's rows, keyed on its header """
        header = self.header
        return [dict(zip(header, line)) for line in self.lines[1:]]

    @property
    def values(self) -> dict:
        """ A key/value section, eg '[Analysis Details]' - first cell to second """
        return {line[0]: (line[1] if len(line) > 1 else None) for line in self.lines if line}


def _clean(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    value = value.strip()
    return value or None


def _split_line(line: str) -> list[str]:
    """ Tab padding is per file rather than per section, so the empty cells off the end of a row
        belong to a wider section and are dropped """
    cells = [c.strip() for c in line.rstrip("\n").split("\t")]
    while cells and not cells[-1]:
        cells.pop()
    return cells


def _is_empty_section(cells: list[str]) -> bool:
    """ 'NA' alone on a row - what a section with nothing in it is written as """
    return cells == [MISSING_VALUE]


def read_combined_variant_output(filename: str) -> dict[str, CombinedVariantOutputSection]:
    """ {section name: section}, in the order the file writes them """
    sections: dict[str, CombinedVariantOutputSection] = {}
    section = None
    with open(filename, encoding="utf-8-sig") as f:
        for line in f:
            cells = _split_line(line)
            if not cells:
                section = None
                continue
            if cells[0].startswith("[") and cells[0].endswith("]"):
                section = CombinedVariantOutputSection(name=cells[0][1:-1])
                sections[section.name] = section
                continue
            if section is None or _is_empty_section(cells):
                continue
            section.lines.append(cells)
    return sections


def can_process_file(filename: str) -> bool:
    """ Recognised from its banner line - a tsv full of gene symbols would otherwise be claimed as a
        gene list """
    try:
        with open(filename, encoding="utf-8-sig") as f:
            first_cells = _split_line(f.readline())
            return bool(first_cells) and bool(FIRST_LINE_PATTERN.fullmatch(first_cells[0]))
    except (OSError, UnicodeDecodeError, IndexError):
        return False


def get_analysis_details(sections: dict[str, CombinedVariantOutputSection]) -> dict:
    """ The pair, its two arms and the software that wrote the file """
    if section := sections.get(ANALYSIS_DETAILS):
        return {k: _clean(v) for k, v in section.values.items()}
    return {}


def get_section_values(sections: dict[str, CombinedVariantOutputSection], name: str) -> dict:
    """ A key/value section, eg '[TMB]' - empty for one this file does not write """
    if section := sections.get(name):
        return {k: _clean(v) for k, v in section.values.items() if k}
    return {}


def get_splice_rows(sections: dict[str, CombinedVariantOutputSection]) -> list[dict]:
    """ The '[Splice Variants]' rows, each exactly as the caller wrote it. A row missing a gene or a
        breakpoint describes nothing we can key on """
    section = sections.get(SPLICE_VARIANTS)
    if section is None:
        return []

    rows = []
    for row in section.rows:
        row = {k: _clean(v) for k, v in row.items() if k}
        if all(row.get(c) for c in SPLICE_COLUMNS):
            rows.append(row)
    return rows


def format_splice_observation(observation: dict) -> str:
    """ One splice row as a reader wants it - the junction and how many reads crossed it """
    parts = []
    breakpoints = [observation.get(BREAKPOINT_1), observation.get(BREAKPOINT_2)]
    if all(breakpoints):
        parts.append("→".join(breakpoints))
    if exon := observation.get(AFFECTED_EXON):
        parts.append(f"exon {exon}")
    if (reads := observation.get(SPLICE_SUPPORTING_READS)) is not None:
        parts.append(f"({reads} reads)")
    return " ".join(parts)
