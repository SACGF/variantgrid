"""
Samtools flagstat output utilities. Flagstat provides simple reads/mapped etc stats on .bam files

@see http://samtools.sourceforge.net/samtools.shtml
"""

from seqauto.qc.qc_utils import match_patterns_in_file


def load_flagstats(flagstats_file):
    PATTERNS = {
        "total": r"^(\d+).*in total",
        "read1": r"^(\d+).*read1",
        "read2": r"^(\d+).*read2",
        "mapped": r"^(\d+).*mapped \(",
        "properly paired": r"^(\d+).*properly paired \(",
    }

    with open(flagstats_file) as f:
        data = match_patterns_in_file(f, PATTERNS, True)
        data = {k: int(match_obj.group(1)) for (k, match_obj) in data.items()}
        return data
