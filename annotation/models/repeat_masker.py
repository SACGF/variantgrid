"""Parse and summarise the VEP RepeatMasker custom annotation.

A structural variant can overlap hundreds of RepeatMasker features. VEP joins the
overlapping subfamily names with ``&`` into a single string (e.g.
``AluSx3&L1M1&MIRb&(CTG)n&...``). Showing that raw is unreadable.

For clinical curation a scientist mostly wants to glance and see *whether* there's repeat
content and roughly what kind, so we roll the subfamily names up to the small set of
RepeatMasker repeat *classes* (SINE/LINE/LTR/...) and count them - the most compressed
view. Classification is name-based (the name is all VEP gives us) so it's best-effort;
anything unmapped falls into ``Other`` (which still signals "repeats are present"), and
the ``&``-separated full list stays available for the exact detail.
"""
import re
from collections import Counter
from dataclasses import dataclass, field

# Repeat class taxonomy and the per-family classifications below follow RepeatMasker /
# Dfam. To verify or extend the mappings, look a family up in Dfam (each family page lists
# its classification, e.g. "LTR/ERV1", "DNA/hAT-Charlie"):
#   - Dfam:        https://www.dfam.org/  (e.g. https://www.dfam.org/family/DF000000001)
#   - RepeatMasker: https://www.repeatmasker.org/
#   - UCSC rmsk repClass field: https://genome.ucsc.edu/cgi-bin/hgTables (group=Repeats)

# RepeatMasker repeat classes (most-compressed grouping for a clinical glance)
SINE = "SINE"
LINE = "LINE"
LTR = "LTR"
DNA = "DNA"
SIMPLE_REPEAT = "Simple repeat"
LOW_COMPLEXITY = "Low complexity"
SATELLITE = "Satellite"
RNA = "RNA"
OTHER = "Other"

# Common human RepeatMasker family prefixes -> class. Kept deliberately small: rare/exotic
# families are left to fall through to OTHER rather than enumerated here.
_CLASS_PREFIXES = (
    (SINE, ("Alu", "MIR", "FLAM", "FRAM", "FAM", "AmnSINE")),
    (LINE, ("L1", "L2", "L3", "L4", "CR1", "HAL1", "RTE")),
    (LTR, ("LTR", "THE1", "MLT", "MST", "HERV", "ERV")),
    (DNA, ("Tigger", "Charlie", "MADE1", "MamRep", "ORSL", "hAT")),
    (RNA, ("5S", "7SK", "7SL", "U6", "tRNA", "rRNA", "snRNA", "scRNA", "srpRNA")),
    (SATELLITE, ("ALR", "HSAT", "GSAT", "SAR", "Satellite", "TAR1")),
)

# MER families split between LTR (ERV) and DNA - the numbering doesn't separate cleanly, so
# the LTR-class numbers are listed and the rest fall through to DNA. Each number below was
# confirmed LTR/ERV via Dfam (see URLs above); unlisted MER numbers default to DNA/hAT.
_LTR_MER_NUMBERS = frozenset({
    "4", "9", "11", "21", "31", "34", "39", "41", "48", "49", "50", "51", "52", "54",
    "57", "61", "65", "66", "67", "68", "70", "72", "73", "74", "76", "77", "83", "84",
    "87", "88", "89", "90", "92", "95", "101", "110",
})
_MER_RE = re.compile(r"^MER(\d+)", re.IGNORECASE)


def classify_repeat(name: str) -> str:
    """ Map a single RepeatMasker subfamily name to its repeat class (best-effort). """
    name = name.strip()
    if not name:
        return OTHER
    if name.startswith("(") and name.endswith(")n"):  # e.g. "(CTG)n", "(AT)n"
        return SIMPLE_REPEAT
    if name.endswith("-rich") or name.lower() in ("polypurine", "polypyrimidine"):
        return LOW_COMPLEXITY
    if m := _MER_RE.match(name):
        return LTR if m.group(1) in _LTR_MER_NUMBERS else DNA
    for repeat_class, prefixes in _CLASS_PREFIXES:
        if name.startswith(prefixes):
            return repeat_class
    if name.endswith("_LINE"):  # e.g. "X7D_LINE"
        return LINE
    return OTHER


@dataclass
class RepeatMaskerSummary:
    """ Grouped summary of a (possibly ``&``-joined) RepeatMasker annotation value. """
    items: list = field(default_factory=list)  # full ordered list, with duplicates
    class_counts: list = field(default_factory=list)  # (class, count) sorted by count desc

    @property
    def total(self) -> int:
        return len(self.items)

    @property
    def is_multiple(self) -> bool:
        return self.total > 1

    @classmethod
    def from_value(cls, value) -> "RepeatMaskerSummary":
        if not value:
            return cls()
        items = [part for part in value.split("&") if part]
        counts = Counter(classify_repeat(name) for name in items)
        # Sort by count descending, then class name for a stable order
        class_counts = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
        return cls(items=items, class_counts=class_counts)
