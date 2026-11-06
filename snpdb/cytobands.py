"""
Cytogenetic band text - a band (11p13), sub-band (11p15.5), arm (17p) or range (11p13-p12, 11p15.5-p15.4) - to
coordinates in a genome build, from the Cytoband rows loaded per contig by snpdb/migrations/0011_cytoband.py.

Entry points: parse_cytoband() (text to CytobandRange, no database) and CytobandRange.get_interval(genome_build).
A band matches its sub-bands by prefix (ISCN numbering is hierarchical: p1 > p13 > p13.1), and a range spans
everything from its first band to its last.
"""
import re
from dataclasses import dataclass
from typing import Optional

from django.db.models import Q

from snpdb.models.models_genome import Contig, Cytoband, GenomeBuild

_CHROM = r"(?:chr)?(?P<{name}>[0-9]{{1,2}}|X|Y)"
_BAND = r"(?P<{name}>[pq](?:\d+(?:\.\d+)?)?)"
CYTOBAND_PATTERN = re.compile(_CHROM.format(name="chrom") + _BAND.format(name="start_band")
                              + r"(?:\s*-\s*(?:" + _CHROM.format(name="end_chrom") + r")?"
                              + _BAND.format(name="end_band") + r")?",
                              re.IGNORECASE)


@dataclass(frozen=True)
class CytobandInterval:
    contig: Contig
    start: int  # 1-based inclusive
    end: int


@dataclass(frozen=True)
class CytobandRange:
    chrom: str
    start_band: str
    end_band: str

    def get_interval(self, genome_build: GenomeBuild) -> Optional[CytobandInterval]:
        """ None unless the chromosome is in the build and both ends of the range are bands on it """
        contig = genome_build.chrom_contig_mappings.get(self.chrom)
        if contig is None:
            return None
        bands = {self.start_band, self.end_band}
        band_q = Q()
        for band in bands:
            band_q |= Q(name__startswith=band)
        cytobands = list(Cytoband.objects.filter(band_q, contig=contig).values_list("name", "start", "end"))
        if not all(any(name.startswith(band) for name, _, _ in cytobands) for band in bands):
            return None
        start = min(start for _, start, _ in cytobands)
        end = max(end for _, _, end in cytobands)
        return CytobandInterval(contig=contig, start=start, end=end)


def parse_cytoband(text: str) -> Optional[CytobandRange]:
    if m := CYTOBAND_PATTERN.fullmatch(text.strip()):
        chrom = m.group("chrom").upper()
        end_chrom = m.group("end_chrom")
        if end_chrom and end_chrom.upper() != chrom:
            return None  # A range across chromosomes isn't one region
        start_band = m.group("start_band").lower()
        end_band = (m.group("end_band") or start_band).lower()
        return CytobandRange(chrom=chrom, start_band=start_band, end_band=end_band)
    return None
