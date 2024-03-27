import re
from dataclasses import dataclass
from typing import Tuple


def format_chrom(chrom, want_chr):
    """ Pass in a chromosome (unknown format), return in your format (CONTIG UPPERCASED)
        @param chrom: chromosome ID (eg 1 or 'chr1')
        @param want_chr: Boolean - whether you want "CHR" at the beginning of chrom
        @return: "chr1" or "1" (for want_chr True/False)
    """

    chrom_no_chr = chrom.upper().replace("CHR", "")
    if want_chr:
        return f"chr{chrom_no_chr}"
    return chrom_no_chr


def get_genomic_size_description(genomic_size):
    if genomic_size == 1000000:
        genomic_size_description = "mb"
    elif genomic_size == 1000:
        genomic_size_description = "kb"
    else:
        genomic_size_description = f"{genomic_size} bases"

    return genomic_size_description


@dataclass
class Range:
    start: int
    end: int


def overlap_fraction(range1: Range, range2: Range) -> float:
    max_start = max(range1.start, range2.start)
    min_end = min(range1.end, range2.end)
    overlap = max(0, min_end - max_start)
    length_range1 = range1.end - range1.start
    if length_range1 == 0:
        return 0
    else:
        return overlap / length_range1


def parse_gnomad_coord(coord: str) -> Tuple[str, int, int]:
    coord_regex = r"(.+):(\d+)-(\d+)"
    m = re.match(coord_regex, coord)  # looks like: chr17:42853981-43253980
    if not m:
        raise ValueError(f"{coord} didn't match regex '{coord_regex}'")
    chrom, start, end = m.groups()
    return chrom, int(start), int(end)
