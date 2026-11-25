"""
Whole-chromosome variants - a chromosome gained or lost (aneuploidy, eg trisomy 21) as one Variant (#2908).

A copy number caller writes a whole-chromosome event as an ordinary SV: DRAGEN's ploidy caller as
`chr21 1 N <DUP> END=46709983` (the contig's length), a segmenting CNV caller as a <DUP> from the first
mappable base to the last. Stored as written these are 46 Mb dups, a different Variant per caller (each
starts and ends somewhere else), and "everyone with trisomy 21" is no one Variant's samples. So the
canonicaliser every insert and lookup goes through (models_variant.py:VariantCoordinate.as_internal_canonical_form)
gives a whole-chromosome gain or loss a single representation:

  * the chromosome's real contig, position 1, ref the fasta base there (N on every standard human contig)
  * alt <CHROM_GAIN> or <CHROM_LOSS> (library/genomics/vcf_enums.py:WHOLE_CHROMOSOME_ALTS) - gain/loss rather
    than trisomy/monosomy because the copy number is the sample's, not the variant's (as GeneLevelSymbolicAlt)
  * svlen = contig length - 1, positive for both, so end = the contig length

So there are exactly two whole-chromosome Variants per chromosome per build. A real contig rather than the
gene-level one because the event has a chromosome: nodes narrow by contig, region filters and IGV read
position and end, and the gene overlaps come from the interval like any other SV too long for VEP.


## What makes a record whole-chromosome

whole_chromosome_alt(): a <DUP> or <DEL> whose span reaches both the first and the last euchromatic band of its
contig - the cytoband file's bands stained gneg or gpos* (models_genome.py:Cytoband.get_euchromatic_bounds). <DUP>
becomes <CHROM_GAIN>, <DEL> <CHROM_LOSS>. <CNV> has no direction and stays as it is, and a contig with no
cytobands (MT, scaffolds, the gene-level contig) has no whole-chromosome form.

Bands rather than a fraction of the contig's length, because a CNV caller starts at the first mappable
sequence: another caller's whole chr21 starts near 13 Mb (the p arm is gvar/stalk/acen bands), only 72% of
the contig, while a 223 Mb chr2 <DUP> ending 18 Mb short (92% of the contig) is not the whole chromosome.
DRAGEN's 1..length records satisfy it trivially. A chromosome arm is still a focal SV.


## What it costs, and where that is paid

The alt is ours, so nothing external can read it: VEP and AnnotSV never see one
(annotation/annotation_pipeline_routing.py, annotation/annotation_version_querysets.py:vep_too_long_q - it gets
a TOO_LONG row with a VariantGeneOverlap per gene and no overlapping_symbols), it has no g.HGVS, c.HGVS or ClinGen
allele (VariantCoordinate.can_be_made_explicit whitelists DEL/DUP/INV), and liftover maps it to the same
chromosome in the other build with no tool (snpdb/liftover.py, AlleleConversionTool.WHOLE_CHROMOSOME).
Variant.get_whole_chromosome_q() / is_whole_chromosome are the predicates.

Searching "trisomy 21", "monosomy X", "chr21 gain" or "loss of chrX" finds the build's whole-chromosome
Variant (parse_whole_chromosome_search, variantopedia/signals/variant_search.py:search_variant_whole_chromosome),
whose samples grid is everyone with that event.
"""
import re
from typing import Optional

from library.genomics.vcf_enums import VCFSymbolicAllele

# The SV a caller writes, to the whole-chromosome alt it becomes
WHOLE_CHROMOSOME_ALT_FOR_SV_ALT = {
    VCFSymbolicAllele.DUP: VCFSymbolicAllele.CHROM_GAIN,
    VCFSymbolicAllele.DEL: VCFSymbolicAllele.CHROM_LOSS,
}

WHOLE_CHROMOSOME_POSITION = 1

_SEARCH_CHROM = r"(?:chr(?:omosome)?\s*)?([0-9]{1,2}|X|Y)"
# "trisomy 21", "monosomy X", "chr21 gain", "gain of chr21"
WHOLE_CHROMOSOME_SEARCH_PATTERN = re.compile(
    rf"^(?:(trisomy|monosomy)\s*{_SEARCH_CHROM}|{_SEARCH_CHROM}\s+(gain|loss)|(gain|loss)\s+of\s+{_SEARCH_CHROM})$",
    re.IGNORECASE)
_SEARCH_WORD_ALTS = {
    "trisomy": VCFSymbolicAllele.CHROM_GAIN,
    "gain": VCFSymbolicAllele.CHROM_GAIN,
    "monosomy": VCFSymbolicAllele.CHROM_LOSS,
    "loss": VCFSymbolicAllele.CHROM_LOSS,
}


def whole_chromosome_alt(alt: str, position: int, end: int,
                         euchromatic_bounds: Optional[tuple[int, int]]) -> Optional[str]:
    """ <CHROM_GAIN>/<CHROM_LOSS> for a <DUP>/<DEL> spanning the contig's euchromatin, else None.
        euchromatic_bounds is Cytoband.get_euchromatic_bounds - None for a contig with no bands """
    whole_alt = WHOLE_CHROMOSOME_ALT_FOR_SV_ALT.get(alt)
    if whole_alt is None or euchromatic_bounds is None:
        return None
    first_band_end, last_band_start = euchromatic_bounds
    if position <= first_band_end and end >= last_band_start:
        return whole_alt
    return None


def whole_chromosome_svlen(contig_length: int) -> int:
    """ Position 1 + svlen is the contig's last base """
    return contig_length - WHOLE_CHROMOSOME_POSITION


def parse_whole_chromosome_search(text: str) -> Optional[tuple[str, str]]:
    """ (chrom, whole-chromosome alt) from a search, eg "trisomy 21" -> ("21", "<CHROM_GAIN>") """
    if m := WHOLE_CHROMOSOME_SEARCH_PATTERN.fullmatch(text.strip()):
        g = m.groups()  # One (word, chrom) pair per alternative in the pattern
        for word, chrom in ((g[0], g[1]), (g[3], g[2]), (g[4], g[5])):
            if word:
                return chrom.upper(), _SEARCH_WORD_ALTS[word.lower()]
    return None

