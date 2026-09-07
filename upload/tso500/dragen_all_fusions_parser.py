"""
Reader for Illumina DRAGEN TSO 500's AllFusions.csv - the rows only, no database work.

What becomes of them is upload.tasks.import_dragen_tso500_all_fusions_task.

The file is a '# key = description' comment block, then a header row, then one row per fusion call
from either of two callers writing into the same file. The header warns their scores and filters are
not comparable, so Caller is carried per row and nothing is keyed on it.

Rows are read unfiltered - 1 of 149 in a real run passes the caller's own filter, and the filter
strings are long semicolon-joined lists, so the same posture as the CNV and small-variant VCFs
applies: take everything, apply our own thresholds later.
"""
import csv
from dataclasses import dataclass, field
from typing import Optional

SOURCE_COMMENT_PREFIX = "# Source = FusionProcessor"
COMMENT_PREFIX = "#"
MISSING_VALUE = "N/A"

CALLER = "Caller"
GENE_A = "Gene A"
GENE_B = "Gene B"
GENE_A_BREAKPOINT = "Gene A Breakpoint"
GENE_B_BREAKPOINT = "Gene B Breakpoint"
DIRECTIONALITY_KNOWN = "Fusion Directionality Known"
ALT_SPLIT = "Alt Split"
ALT_PAIR = "Alt Pair"
ALT_SPLIT_DEDUP = "Alt Split Dedup"
ALT_PAIR_DEDUP = "Alt Pair Dedup"

REQUIRED_COLUMNS = (CALLER, GENE_A, GENE_B, GENE_A_BREAKPOINT, GENE_B_BREAKPOINT)

# INFO fields the rows are carried in once they become a VCF - they land in CohortGenotype.info via
# the standard bulk importer, which stores every INFO field the header declares
FUSION_INFO = "FUSION"
FUSION_OBSERVATIONS_INFO = "FUSION_OBS"
# What separates one call from the next in the text the grid and the merge record show. The parts of
# a call - caller, breakpoints, read count - never contain it
OBSERVATION_SEPARATOR = "; "


@dataclass
class AllFusionsRow:
    """ One call. `data` keeps every column as written, which is what reaches CohortGenotype.info """
    caller: str
    gene_a: str
    gene_b: str
    data: dict = field(default_factory=dict)

    @property
    def gene_a_breakpoint(self) -> Optional[str]:
        return self.data.get(GENE_A_BREAKPOINT)

    @property
    def gene_b_breakpoint(self) -> Optional[str]:
        return self.data.get(GENE_B_BREAKPOINT)

    @property
    def directionality_known(self) -> bool:
        """ Reported by the caller rather than inferred from gene order """
        return str(self.data.get(DIRECTIONALITY_KNOWN, "")).strip().lower() == "true"


def _clean(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    value = value.strip()
    if not value or value == MISSING_VALUE:
        return None
    return value


def _read_header_and_rows(f) -> tuple[list[str], csv.DictReader]:
    comments = []
    for line in f:
        if not line.startswith(COMMENT_PREFIX):
            # csv.reader over the remaining handle, with this line put back as the header
            reader = csv.DictReader([line] + f.readlines())
            return comments, reader
        comments.append(line.rstrip("\n"))
    return comments, csv.DictReader([])


def can_process_file(filename: str) -> bool:
    """ Recognised from its own '# Source = FusionProcessor' line plus the header row - a gene list
        would otherwise claim it, being a csv full of gene symbols """
    try:
        with open(filename, "rt", encoding="utf-8-sig") as f:
            first_line = f.readline()
            if not first_line.startswith(SOURCE_COMMENT_PREFIX):
                return False
            _comments, reader = _read_header_and_rows(f)
            fieldnames = reader.fieldnames or []
            return all(c in fieldnames for c in REQUIRED_COLUMNS)
    except (OSError, UnicodeDecodeError, csv.Error):
        return False


def read_all_fusions(filename: str) -> tuple[list[str], list[AllFusionsRow]]:
    """ (the file's comment block, one AllFusionsRow per call) """
    with open(filename, "rt", encoding="utf-8-sig") as f:
        comments, reader = _read_header_and_rows(f)
        if missing := [c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]:
            raise ValueError(f"{filename} is missing column(s): {', '.join(missing)}")

        rows = []
        for record in reader:
            data = {k: _clean(v) for k, v in record.items() if k is not None}
            gene_a = data.get(GENE_A)
            gene_b = data.get(GENE_B)
            if gene_a is None and gene_b is None:
                continue  # Named no genes at all, so describes nothing we can key on
            rows.append(AllFusionsRow(caller=data.get(CALLER), gene_a=gene_a, gene_b=gene_b, data=data))
    return comments, rows


def _supporting_reads(observation: dict) -> Optional[int]:
    """ Reads supporting the fusion - the caller's deduplicated counts where it reported them """
    for split_key, pair_key in ((ALT_SPLIT_DEDUP, ALT_PAIR_DEDUP), (ALT_SPLIT, ALT_PAIR)):
        values = [observation.get(split_key), observation.get(pair_key)]
        if any(v is not None for v in values):
            try:
                return sum(int(float(v)) for v in values if v is not None)
            except ValueError:
                return None
    return None


def format_fusion_observation(observation: dict) -> str:
    """ One caller row as a reader wants it - who called it, the breakpoints, how many reads """
    parts = [observation.get(CALLER) or "?"]
    breakpoints = [observation.get(GENE_A_BREAKPOINT), observation.get(GENE_B_BREAKPOINT)]
    if all(breakpoints):
        parts.append("\u2192".join(breakpoints))
    if (reads := _supporting_reads(observation)) is not None:
        parts.append(f"({reads} reads)")
    return " ".join(parts)


def format_fusion_observations(observations) -> str:
    """ Every call a gene pair was merged from, in one line """
    return OBSERVATION_SEPARATOR.join(format_fusion_observation(o) for o in observations)
