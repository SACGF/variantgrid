"""
Reader for Illumina DRAGEN TSO 500's MetricsOutput.tsv - the file only, no database work.

What becomes of it is upload.tasks.import_dragen_tso500_metrics_output_task, via
upload.tso500.dragen_metrics_output_records. The file is written once per run beside the
CombinedVariantOutputs and is that file's layout - a banner, then one '[Section]' of key/values or a
table - so the CVO's section reader reads it whole (read_sections).

Entry points are can_process_file (the banner line, since a tsv would otherwise be claimed as a gene
list), read_sections and read_library_qc, which turns the QC sections into one LibraryQCMetrics per
column.

A column is a pair, named by the CVO's 'Pair ID', and carries both of its arms: the DNA sections and
the RNA sections have values in the same column, with 'NA' down an arm the pair does not have. A
metric passes when LSL <= value <= USL, an 'NA' guideline being no bound - the same rule the lab's
own release summary applies. The file's own guideline is the policy, so there are no threshold
settings; settings.TSO500_LIBRARY_QC_GUIDELINES is the exception, for the one metric a lab judges by
its own number, and each metric records which of the two it was judged by (guideline_source).

Section names and metric lists move between module versions (2.6.2 adds PCT_CHIMERIC_READS to the
small-variant section and EXCESSIVE_TF to GIS, and drops PCT_PF_UQ_READS), so only the section names
are known here: every metric of a known section counts towards its category, and a section this has
no category for is skipped. '[Run QC Metrics]' describes the run rather than any library, and the
'[... Expanded Metrics]' sections carry no guidelines, so neither is read.
"""
import logging
import re
from dataclasses import dataclass, field
from typing import Optional

from django.conf import settings

from patients.models_enums import NucleicAcid
from seqauto.models.models_enums import LibraryQCCategory
from upload.tso500.dragen_combined_variant_output_parser import (
    MISSING_VALUE,
    CombinedVariantOutputSection,
    read_combined_variant_output,
)

# The format is the CVO's, so the same reader does both - named for what it does here
read_sections = read_combined_variant_output

# The banner line, which is what recognises the file. 2.6 puts its version in it, the way the CVO's does
FIRST_LINE_PATTERN = re.compile(
    r"DRAGEN TruSight Oncology 500 (?:v[\d.]+ )?Analysis Software - Metrics Output")

HEADER = "Header"
ANALYSIS_STATUS = "Analysis Status"

# [Header] keys - Output Date / Output Time are the CVO's, so measured_date reads both
WORKFLOW_VERSION = "Workflow Version"

# [Analysis Status] rows. The header's first cell is empty, so the sample names start at column two
COMPLETED_ALL_STEPS = "COMPLETED_ALL_STEPS"

# The QC sections' first three columns - everything after them is a sample
METRIC = "Metric (UOM)"
LSL_GUIDELINE = "LSL Guideline"
USL_GUIDELINE = "USL Guideline"
GUIDELINE_COLUMNS = 3

METRIC_UNIT_PATTERN = re.compile(r"^(?P<metric>\S+)\s*\((?P<unit>[^)]*)\)$")

# Which library QC category each section vouches for. What DRAGEN calls them is this module's
# business; what they mean is LibraryQCCategory's
SECTION_CATEGORIES = {
    "DNA Library QC Metrics": LibraryQCCategory.DNA,
    "DNA Library QC Metrics for Small Variant Calling and TMB": LibraryQCCategory.SMALL_VARIANT_TMB,
    "DNA Library QC Metrics for MSI": LibraryQCCategory.MSI,
    "DNA Library QC Metrics for CNV": LibraryQCCategory.CNV,
    "DNA Library QC Metrics for GIS": LibraryQCCategory.GIS,
    "RNA Library QC Metrics": LibraryQCCategory.RNA,
}

# Which arm each category describes - the file's column is the pair, so the row says which half of it
CATEGORY_NUCLEIC_ACID = {
    LibraryQCCategory.DNA: NucleicAcid.DNA,
    LibraryQCCategory.SMALL_VARIANT_TMB: NucleicAcid.DNA,
    LibraryQCCategory.MSI: NucleicAcid.DNA,
    LibraryQCCategory.CNV: NucleicAcid.DNA,
    LibraryQCCategory.GIS: NucleicAcid.DNA,
    LibraryQCCategory.RNA: NucleicAcid.RNA,
}

GUIDELINE_SETTING = "TSO500_LIBRARY_QC_GUIDELINES"
FILE_GUIDELINE_SOURCE = "file"


@dataclass
class LibraryQCMetrics:
    """ One column of the file - what every QC section said about that pair, both arms """
    pair_id: str
    completed: Optional[bool] = None
    # {LibraryQCCategory: {metric: {"value", "unit", "lsl", "usl", "passed"}}}
    categories: dict[str, dict] = field(default_factory=dict)

    def category_passed(self, category: str) -> Optional[bool]:
        """ All of the category's metrics were within their guideline - None where every value in
            the section is NA, which is what the arm this library is not looks like """
        outcomes = [metric["passed"] for metric in self.categories.get(category, {}).values()]
        outcomes = [passed for passed in outcomes if passed is not None]
        if not outcomes:
            return None
        return all(outcomes)


def can_process_file(filename: str) -> bool:
    """ Recognised from its banner line - a tsv of gene symbols would otherwise be claimed as a gene list """
    try:
        with open(filename, encoding="utf-8-sig") as f:
            first_cell = f.readline().split("\t")[0].strip()
            return bool(FIRST_LINE_PATTERN.fullmatch(first_cell))
    except (OSError, UnicodeDecodeError, IndexError):
        return False


def _number(value: Optional[str]) -> Optional[float]:
    """ A cell as a number - 'NA' (no bound, or an arm that was not sequenced) is None """
    if value is None or value.strip() in ("", MISSING_VALUE):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_metric(cell: str) -> tuple[str, Optional[str]]:
    """ 'PCT_EXON_50X (%)' as its name and its unit - '(NA)' is a metric with no unit """
    if m := METRIC_UNIT_PATTERN.match(cell):
        unit = m.group("unit")
        return m.group("metric"), None if unit == MISSING_VALUE else unit
    return cell, None


def metric_passed(value: Optional[float], lsl: Optional[float], usl: Optional[float]) -> Optional[bool]:
    """ Within its guideline, an absent bound being no bound. No value is no judgement """
    if value is None:
        return None
    if lsl is not None and value < lsl:
        return False
    if usl is not None and value > usl:
        return False
    return True


def get_workflow_version(sections: dict[str, CombinedVariantOutputSection]) -> Optional[str]:
    if section := sections.get(HEADER):
        return section.values.get(WORKFLOW_VERSION)
    return None


def _pair_columns(header: list[str]) -> list[str]:
    return [name for name in header[GUIDELINE_COLUMNS:] if name]


def _analysis_status(sections: dict[str, CombinedVariantOutputSection]) -> dict[str, Optional[bool]]:
    """ COMPLETED_ALL_STEPS per pair. The header's first cell is empty, so the names start at
        column two - a pair whose run did not finish is what the report's Fail caveat is about """
    section = sections.get(ANALYSIS_STATUS)
    if section is None or not section.lines:
        return {}
    pair_ids = section.lines[0][1:]
    completed = {}
    for line in section.lines[1:]:
        if line[0] != COMPLETED_ALL_STEPS:
            continue
        for pair_id, value in zip(pair_ids, line[1:]):
            if pair_id:
                completed[pair_id] = {"TRUE": True, "FALSE": False}.get(value.upper())
    return completed


def guideline_for(section_name: str, metric: str, lsl: Optional[float],
                  usl: Optional[float]) -> tuple[Optional[float], Optional[float], str]:
    """ The bounds a metric is judged by and where they came from. The file's own guideline is the
        policy unless the lab has set its own for that metric - keyed on (section, metric), since
        MEDIAN_INSERT_SIZE appears in two sections with different numbers """
    overrides = getattr(settings, GUIDELINE_SETTING, None) or {}
    if override := overrides.get((section_name, metric)):
        return override[0], override[1], GUIDELINE_SETTING
    return lsl, usl, FILE_GUIDELINE_SOURCE


def read_library_qc(sections: dict[str, CombinedVariantOutputSection]) -> list[LibraryQCMetrics]:
    """ One LibraryQCMetrics per column, in the order the file writes them """
    completed = _analysis_status(sections)
    libraries: dict[str, LibraryQCMetrics] = {
        pair_id: LibraryQCMetrics(pair_id=pair_id, completed=is_completed)
        for pair_id, is_completed in completed.items()
    }

    for section_name, category in SECTION_CATEGORIES.items():
        section = sections.get(section_name)
        if section is None:
            logging.info("MetricsOutput has no '[%s]' section", section_name)
            continue
        pair_ids = _pair_columns(section.header)
        for line in section.lines[1:]:
            metric, unit = parse_metric(line[0])
            lsl, usl, guideline_source = guideline_for(
                section_name, metric,
                _number(line[1] if len(line) > 1 else None),
                _number(line[2] if len(line) > 2 else None))
            for pair_id, cell in zip(pair_ids, line[GUIDELINE_COLUMNS:]):
                value = _number(cell)
                library = libraries.setdefault(pair_id, LibraryQCMetrics(pair_id=pair_id))
                library.categories.setdefault(category, {})[metric] = {
                    "value": value,
                    "unit": unit,
                    "lsl": lsl,
                    "usl": usl,
                    "passed": metric_passed(value, lsl, usl),
                    "guideline_source": guideline_source,
                }

    return list(libraries.values())
