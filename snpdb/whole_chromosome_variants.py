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
  * svlen = contig length - 1, positive for all three alts, so end = the contig length

So there are exactly three whole-chromosome Variants per chromosome per build, the third the euploid one below. A real contig rather than the
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

## Euploid chromosomes (#2915)

A ploidy caller writes each euploid chromosome as a reference span (ALT '.', END the contig's end, no GT).
That is a measurement - the chromosome is at the expected copy number for the sample's karyotype, so a male's
single chrX is euploid - and "tested and normal for chr21" is worth keeping. A reference span reaching the same
slack (is_whole_chromosome_span) is stored as <CNV:CHROM>, the spec's CNV being a "region of uniform copy
number (may be deletion, duplication or copy number neutral)": vcf_clean_and_filter keeps it where it skips
any other reference span, and the importer's whole_chromosome_reference_alt reads it as <CNV:CHROM>. Its
genotype is HOM_REF, as for any ALT '.' record in a VCF without GT
(upload/vcf/bulk_genotype_vcf_processor.py), so SampleNode shows it only with ref ticked
(analysis/models/nodes/sources/sample_node.py:SampleNode._get_zygosities) and it is not a call in the sample's stats.

A caller's whole-chromosome <CNV> has no direction and is stored as written, so <CNV:CHROM> is only ever a
reference call.

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

It is named by ISCN's numerical abnormality, the chromosome with + or - (whole_chromosome_iscn: "+21", "-X"),
as HGVS has no form for it (its complex-variant recommendations defer chromosomal changes to ISCN). With no copy
number or sex this is the variant-level form: a constitutional monosomy X is "45,X" in a karyotype, which is the
sample's. Grids, the variant page title and its card show whole_chromosome_label, "+21 GAIN" (#2910).

Searching "trisomy 21", "monosomy X", "chr21 gain", "loss of chrX" or the ISCN ("+21", or as labelled "+21 GAIN")
finds the build's whole-chromosome Variant (parse_whole_chromosome_search,
variantopedia/signals/variant_search.py:search_variant_whole_chromosome), whose samples grid is everyone with
that event.
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

# ISCN's sign for a whole chromosome gained or lost, and the word the label and badge say. A euploid
# chromosome has no ISCN abnormality, so is labelled by its contig
_ISCN_SIGNS = {
    VCFSymbolicAllele.DUP_CHROM: "+",
    VCFSymbolicAllele.DEL_CHROM: "-",
}
_LABEL_WORDS = {
    VCFSymbolicAllele.DUP_CHROM: "GAIN",
    VCFSymbolicAllele.DEL_CHROM: "LOSS",
    VCFSymbolicAllele.CNV_CHROM: "NORMAL",
}

_SEARCH_CHROM = r"(?:chr(?:omosome)?\s*)?([0-9]{1,2}|X|Y)"
# "trisomy 21", "monosomy X", "chr21 gain", "gain of chr21", the ISCN "+21" / "-X" alone or as labelled
# ("+21 GAIN"), and a euploid one as labelled ("chr21 NORMAL")
WHOLE_CHROMOSOME_SEARCH_PATTERN = re.compile(
    rf"^(?:(trisomy|monosomy)\s*{_SEARCH_CHROM}|{_SEARCH_CHROM}\s+(gain|loss|normal)|(gain|loss)\s+of\s+{_SEARCH_CHROM}"
    rf"|([+-])\s*([0-9]{{1,2}}|X|Y)(?:\s+(?:gain|loss))?)$",
    re.IGNORECASE)
_SEARCH_WORD_ALTS = {
    "trisomy": VCFSymbolicAllele.DUP_CHROM,
    "gain": VCFSymbolicAllele.DUP_CHROM,
    "+": VCFSymbolicAllele.DUP_CHROM,
    "monosomy": VCFSymbolicAllele.DEL_CHROM,
    "loss": VCFSymbolicAllele.DEL_CHROM,
    "-": VCFSymbolicAllele.DEL_CHROM,
    "normal": VCFSymbolicAllele.CNV_CHROM,
}


def _uncalled_ends(position: int, end: int, span: tuple[int, int]) -> tuple[int, int]:
    """ How far short of each end of span a record stops """
    span_start, span_end = span
    return max(0, position - span_start), max(0, span_end - end)


def is_whole_chromosome_span(position: int, end: int, span: Optional[tuple[int, int]]) -> bool:
    """ position..end reaches within VARIANT_WHOLE_CHROMOSOME_END_SLACK of each end of span.
        span is Cytoband.get_whole_chromosome_span - None for a contig with no bands """
    if span is None:
        return False
    slack = settings.VARIANT_WHOLE_CHROMOSOME_END_SLACK
    return all(uncalled <= slack for uncalled in _uncalled_ends(position, end, span))


def whole_chromosome_alt(alt: str, position: int, end: int, span: Optional[tuple[int, int]]) -> Optional[str]:
    """ <DUP:CHROM>/<DEL:CHROM> for a <DUP>/<DEL> spanning the whole chromosome, else None """
    whole_alt = WHOLE_CHROMOSOME_ALT_FOR_SV_ALT.get(alt)
    if whole_alt and is_whole_chromosome_span(position, end, span):
        return whole_alt
    return None


def whole_chromosome_reference_alt(position: int, end: int, span: Optional[tuple[int, int]]) -> Optional[str]:
    """ <CNV:CHROM> for a reference span (ALT '.' with END) over the whole chromosome - a euploid one - else None """
    if is_whole_chromosome_span(position, end, span):
        return VCFSymbolicAllele.CNV_CHROM
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


def whole_chromosome_iscn(chrom: str, alt: str) -> Optional[str]:
    """ ISCN for a whole-chromosome alt, eg ("21", "<DUP:CHROM>") -> "+21", ("chrX", "<DEL:CHROM>") -> "-X".
        None for any other alt """
    if sign := _ISCN_SIGNS.get(alt):
        return sign + re.sub(r"^chr", "", chrom, flags=re.IGNORECASE)
    return None


def whole_chromosome_label(chrom: str, alt: str) -> Optional[str]:
    """ "+21 GAIN", "chr21 NORMAL" - as the grid's label and badge (variantgrid_formats.js:WHOLE_CHROMOSOME_KINDS) """
    if not (word := _LABEL_WORDS.get(alt)):
        return None
    name = whole_chromosome_iscn(chrom, alt) or ("chr" + re.sub(r"^chr", "", chrom, flags=re.IGNORECASE))
    return f"{name} {word}"


def parse_whole_chromosome_search(text: str) -> Optional[tuple[str, str]]:
    """ (chrom, whole-chromosome alt) from a search, eg "trisomy 21" -> ("21", "<DUP:CHROM>") """
    if m := WHOLE_CHROMOSOME_SEARCH_PATTERN.fullmatch(text.strip()):
        g = m.groups()  # One (word, chrom) pair per alternative in the pattern
        for word, chrom in ((g[0], g[1]), (g[3], g[2]), (g[4], g[5]), (g[6], g[7])):
            if word:
                return chrom.upper(), _SEARCH_WORD_ALTS[word.lower()]
    return None

