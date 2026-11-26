"""
Whole-chromosome variants - a chromosome gained or lost (aneuploidy, eg trisomy 21) as one Variant (#2908).

A copy number caller writes a whole-chromosome event as an ordinary SV: DRAGEN's ploidy caller as
`chr21 1 N <DUP> END=46709983` (the contig's length), a segmenting CNV caller as a <DUP> from the first
mappable base to the last. Stored as written these are 46 Mb dups, a different Variant per caller (each
starts and ends somewhere else), and "everyone with trisomy 21" is no one Variant's samples. So the
canonicaliser every insert and lookup goes through (models_variant.py:VariantCoordinate.as_internal_canonical_form)
gives a whole-chromosome gain or loss a single representation:

  * the chromosome's real contig, position 1, ref the fasta base there (N on every standard human contig)
  * alt <DUP:CHROM> or <DEL:CHROM> (library/genomics/vcf_enums.py:WHOLE_CHROMOSOME_ALTS) - gain/loss rather
    than trisomy/monosomy because the copy number is the sample's, not the variant's (as GeneLevelSymbolicAlt)
  * svlen = contig length - 1, positive for both, so end = the contig length

So there are exactly two whole-chromosome Variants per chromosome per build. A real contig rather than the
gene-level one because the event has a chromosome: nodes narrow by contig, region filters and IGV read
position and end, and the gene overlaps come from the interval like any other SV too long for VEP.


## What makes a record whole-chromosome

whole_chromosome_alt(): a <DUP> or <DEL> reaching within settings.VARIANT_WHOLE_CHROMOSOME_END_SLACK (1 Mb) of
each end of the chromosome's whole-chromosome span (models_genome.py:Cytoband.get_whole_chromosome_span): from the
start of its first euchromatic band (stained gneg or gpos*) to the end of its last, less chrY's pseudoautosomal
regions. <DUP> becomes <DUP:CHROM>, <DEL> <DEL:CHROM>. <CNV> has no direction and stays as it is, and a contig with
no cytobands (MT, scaffolds, the gene-level contig) has no whole-chromosome form.

The span is what callers cover. DRAGEN's ploidy caller writes 1..contig length. A segmenting caller starts at its
first usable bin: past the telomere's N on most chromosomes, anywhere from the first non-N base to the start of
the q arm on an acrocentric (chr21 has sequence from 5 Mb, euchromatin from 13 Mb), and outside the PARs on chrY,
which analysis-set references mask. The slack is smaller than the terminal bands, so a call missing one
(eg a chr21 dup ending before q22.3) is a partial gain or loss, as is any arm.

A <DUP>/<DEL> that misses the slack but covers VARIANT_WHOLE_CHROMOSOME_NEAR_MISS_FRACTION (90%) of the span
(is_whole_chromosome_near_miss) is stored as written, and the VCF import says it may have been a whole
chromosome (upload/vcf/abstract_bulk_vcf_processor.py) and reports it to Rollbar. If it was, raising the slack
in the deployment's settings and running one_off_whole_chromosome_variants moves it.


## What it costs, and where that is paid

The alts are our own subtypes of the VCF spec's DUP/DEL ("Implementations are free to define their own
subtypes"), so a reader of an exported VCF sees a duplication or deletion. Nothing external annotates them: VEP
and AnnotSV never see one
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

from django.conf import settings

from library.genomics.vcf_enums import VCFSymbolicAllele

# The SV a caller writes, to the whole-chromosome alt it becomes
WHOLE_CHROMOSOME_ALT_FOR_SV_ALT = {
    VCFSymbolicAllele.DUP: VCFSymbolicAllele.DUP_CHROM,
    VCFSymbolicAllele.DEL: VCFSymbolicAllele.DEL_CHROM,
}

WHOLE_CHROMOSOME_POSITION = 1

_SEARCH_CHROM = r"(?:chr(?:omosome)?\s*)?([0-9]{1,2}|X|Y)"
# "trisomy 21", "monosomy X", "chr21 gain", "gain of chr21"
WHOLE_CHROMOSOME_SEARCH_PATTERN = re.compile(
    rf"^(?:(trisomy|monosomy)\s*{_SEARCH_CHROM}|{_SEARCH_CHROM}\s+(gain|loss)|(gain|loss)\s+of\s+{_SEARCH_CHROM})$",
    re.IGNORECASE)
_SEARCH_WORD_ALTS = {
    "trisomy": VCFSymbolicAllele.DUP_CHROM,
    "gain": VCFSymbolicAllele.DUP_CHROM,
    "monosomy": VCFSymbolicAllele.DEL_CHROM,
    "loss": VCFSymbolicAllele.DEL_CHROM,
}


def _uncalled_ends(position: int, end: int, span: tuple[int, int]) -> tuple[int, int]:
    """ How far short of each end of span a record stops """
    span_start, span_end = span
    return max(0, position - span_start), max(0, span_end - end)


def whole_chromosome_alt(alt: str, position: int, end: int, span: Optional[tuple[int, int]]) -> Optional[str]:
    """ <DUP:CHROM>/<DEL:CHROM> for a <DUP>/<DEL> reaching within VARIANT_WHOLE_CHROMOSOME_END_SLACK of each end
        of span, else None. span is Cytoband.get_whole_chromosome_span - None for a contig with no bands """
    whole_alt = WHOLE_CHROMOSOME_ALT_FOR_SV_ALT.get(alt)
    if whole_alt is None or span is None:
        return None
    slack = settings.VARIANT_WHOLE_CHROMOSOME_END_SLACK
    if all(uncalled <= slack for uncalled in _uncalled_ends(position, end, span)):
        return whole_alt
    return None


def is_whole_chromosome_near_miss(alt: str, position: int, end: int, span: Optional[tuple[int, int]]) -> bool:
    """ A <DUP>/<DEL> that isn't whole-chromosome but covers VARIANT_WHOLE_CHROMOSOME_NEAR_MISS_FRACTION of span -
        possibly a caller's whole chromosome that the end slack is too tight for """
    if alt not in WHOLE_CHROMOSOME_ALT_FOR_SV_ALT or span is None or whole_chromosome_alt(alt, position, end, span):
        return False
    span_start, span_end = span
    span_length = span_end - span_start
    covered = span_length - sum(_uncalled_ends(position, end, span))
    return covered >= settings.VARIANT_WHOLE_CHROMOSOME_NEAR_MISS_FRACTION * span_length


def whole_chromosome_svlen(contig_length: int) -> int:
    """ Position 1 + svlen is the contig's last base """
    return contig_length - WHOLE_CHROMOSOME_POSITION


def parse_whole_chromosome_search(text: str) -> Optional[tuple[str, str]]:
    """ (chrom, whole-chromosome alt) from a search, eg "trisomy 21" -> ("21", "<DUP:CHROM>") """
    if m := WHOLE_CHROMOSOME_SEARCH_PATTERN.fullmatch(text.strip()):
        g = m.groups()  # One (word, chrom) pair per alternative in the pattern
        for word, chrom in ((g[0], g[1]), (g[3], g[2]), (g[4], g[5])):
            if word:
                return chrom.upper(), _SEARCH_WORD_ALTS[word.lower()]
    return None

