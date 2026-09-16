"""
Splice calls - "AR-V7", "MET exon 14 skipping" - as gene-level variants.

@see snpdb.gene_level_variants for why one of these is a Variant, and genes.gene_level_resolver for
turning the caller's gene name into the identity it is stored under.

Identity is the gene plus the junction's label, so two events in one gene are two variants. The
label comes from genes.models.models_splice_event.SpliceEvent, which names the junctions a report
talks about; a junction we have no name for is labelled with its own coordinates
(X_66905968_66914514), which reads as raw coordinates on a report and is the prompt to add a row.

A splice event has no record of its own the way a fusion has a GeneFusion - the alt carries the gene
and the label, and the caller's breakpoints ride along in the VCF's INFO - so the objects here are
resolved identities and a read-only view of a Variant, not models.
"""
import re
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from typing import Optional

from genes.gene_level_resolver import GeneLevelNameResolver
from genes.models import GeneLevelId
from genes.models.models_splice_event import SpliceEvent
from library.genomics.vcf_enums import GeneLevelSymbolicAlt
from snpdb.gene_level_variants import (
    GENE_LEVEL_CONTIG_NAME,
    GENE_LEVEL_REF,
    GENE_LEVEL_SVLEN,
)
from snpdb.models import Contig, GenomeBuild, Variant, VariantCoordinate

# 'X_66905968_66914514' - the label a junction we have no SpliceEvent for is stored under. The
# contig is spelled out rather than left as \w+, so 'ARX_66905968_66914514' splits after the gene
COORDINATE_LABEL_PATTERN = r"(?:chr)?(?:[0-9]{1,2}|MT|X|Y|M)_[0-9]+_[0-9]+"
# The label shapes a junction is named with, after the gene: 'V7', 'vIII', 'ex14skip',
# 'exon 14 skipping', or the coordinate label. Widen this when a new label shape is seeded, or
# search will not offer it - resolution itself matches on the name, not on this.
SPLICE_LABEL_PATTERN = rf"v[0-9]+|v[IVX]+|ex(?:on)?\s*[0-9]+\s*skip(?:ping)?|{COORDINATE_LABEL_PATTERN}"
# 'AR V7', 'AR-V7 splice variant', 'MET exon 14 skipping'. The space between gene and label is
# optional because a classification's imported c.HGVS reaches us with its spaces removed
# (@see ImportedAlleleInfo._tidy_input_value); a gene that has to resolve is what keeps this from
# claiming ordinary words.
SPLICE_STRING_PATTERN = re.compile(
    rf"^\s*[A-Za-z0-9.]+?\s*[-_]?\s*(?:{SPLICE_LABEL_PATTERN})(?:\s*splice(?:\s*variant)?)?\s*$",
    re.IGNORECASE)
# The coordinate form is the one written form we split, since no row names its junction
SPLICE_COORDINATE_PATTERN = re.compile(
    rf"^\s*(?P<gene>[A-Za-z0-9.]+?)\s*[-_]?\s*(?P<label>{COORDINATE_LABEL_PATTERN})\s*$", re.IGNORECASE)
# Everything a written form may put between its words, so one string is one key
SPLICE_KEY_SEPARATORS = re.compile(r"[\s\-_]")


def splice_canonical_str(gene: GeneLevelId, label: str) -> str:
    """ 'AR V7' - the gene then the junction's label, the form to display and to send anywhere off
        this deployment (@see GeneLevelId). SpliceEvent.display is what a report writes instead,
        where we have a row for the junction. """
    return f"{gene.symbol_str} {label}"


def coordinate_label(contig: Contig, donor: int, acceptor: int) -> str:
    """ The label for a junction no SpliceEvent names. Build-specific by construction, and it reads
        as raw coordinates wherever the label is printed - which is the prompt to add a row """
    return f"{contig.name.replace('chr', '')}_{donor}_{acceptor}"


@dataclass(frozen=True)
class ResolvedSpliceEvent:
    """ A splice event's identity, before it has a Variant. The gene is the Locus.position and the
        alt carries it again plus the label - @see snpdb.gene_level_variants """
    gene: GeneLevelId
    label: str
    splice_event: Optional[SpliceEvent] = None  # The junction's name, where we have one

    @property
    def alt(self) -> str:
        return GeneLevelSymbolicAlt.format(GeneLevelSymbolicAlt.SPLICE, self.gene.alt_namespace,
                                           self.gene.pk, self.label)

    @property
    def variant_coordinate(self) -> VariantCoordinate:
        """ What the loader writes as a VCF record, and looks the Variant up by afterwards """
        return VariantCoordinate(chrom=GENE_LEVEL_CONTIG_NAME, position=self.gene.pk,
                                 ref=GENE_LEVEL_REF, alt=self.alt, svlen=GENE_LEVEL_SVLEN)

    @property
    def canonical_str(self) -> str:
        return splice_canonical_str(self.gene, self.label)

    @property
    def display(self) -> str:
        """ What a report writes - the named junction where we have it, else the canonical form """
        if self.splice_event:
            return self.splice_event.display
        return self.canonical_str


class SpliceEventResolver:
    """ Holds the symbol caches and the build's contigs, so build one per file rather than one per row """

    def __init__(self, genome_build: GenomeBuild, name_resolver: GeneLevelNameResolver = None):
        self.genome_build = genome_build
        self.name_resolver = name_resolver or GeneLevelNameResolver()

    def get_splice_event(self, contig: Contig, donor: int, acceptor: int) -> Optional[SpliceEvent]:
        return SpliceEvent.objects.filter(genome_build=self.genome_build, contig=contig,
                                          donor=donor, acceptor=acceptor).first()

    def resolve(self, gene_name: str, chrom: str, donor: int, acceptor: int) -> Optional[ResolvedSpliceEvent]:
        """ The identity a caller's row is stored under. The junction's coordinates name it where we
            have a SpliceEvent for them; otherwise the coordinates themselves are the label """

        resolved_gene = self.name_resolver.resolve_gene(gene_name)
        if resolved_gene is None:
            return None

        contig = self.genome_build.chrom_contig_mappings.get(chrom)
        if contig is None:
            return None

        splice_event = self.get_splice_event(contig, donor, acceptor)
        label = splice_event.label if splice_event else coordinate_label(contig, donor, acceptor)
        return ResolvedSpliceEvent(gene=resolved_gene.gene_level_id, label=label,
                                   splice_event=splice_event)


@dataclass(frozen=True)
class SpliceEventVariant:
    """ A splice Variant read back off its alt - what a fusion gets from its GeneFusion row. The alt
        holds the whole identity, so nothing else has to be stored or joined """
    variant: Variant
    gene: GeneLevelId
    label: str

    @property
    def gene_level_ids(self) -> list[GeneLevelId]:
        """ The genes this event is about, so gene lists and annotation read it the way they read a
            fusion's partners """
        return [self.gene]

    @property
    def canonical_str(self) -> str:
        return splice_canonical_str(self.gene, self.label)

    @property
    def splice_event(self) -> Optional[SpliceEvent]:
        """ The junction's name. Keyed on (gene symbol, label) rather than the coordinates, which
            the Variant does not carry - the pair is unique per build, and a label means the same
            event in every build """
        return SpliceEvent.objects.filter(gene_symbol=self.gene.gene_symbol_id,
                                          label=self.label).first()

    @property
    def display(self) -> str:
        if splice_event := self.splice_event:
            return splice_event.display
        return self.canonical_str


def get_splice_event_variant(variant: Variant) -> Optional[SpliceEventVariant]:
    """ The splice event a gene-level Variant is, or None if it is a fusion / copy number call """
    if not variant.is_gene_level:
        return None
    parsed = GeneLevelSymbolicAlt.parse(variant.alt.seq)
    if parsed is None:
        return None
    kind, _namespace, gene_id, label = parsed
    if kind != GeneLevelSymbolicAlt.SPLICE:
        return None
    gene = GeneLevelId.objects.filter(pk=gene_id).first()
    if gene is None:
        return None
    return SpliceEventVariant(variant=variant, gene=gene, label=label)


def splice_event_variants(variant_qs) -> Iterator[SpliceEventVariant]:
    """ Every splice event among gene-level variants - what annotation walks, in place of the rows a
        fusion or a copy number call has """
    variants: Iterable[Variant] = variant_qs.filter(Variant.get_gene_level_q()) \
                                            .select_related("locus", "alt").iterator()
    gene_ids: dict[int, GeneLevelId] = {}
    for variant in variants:
        parsed = GeneLevelSymbolicAlt.parse(variant.alt.seq)
        if parsed is None:
            continue
        kind, _namespace, gene_id, label = parsed
        if kind != GeneLevelSymbolicAlt.SPLICE:
            continue
        gene = gene_ids.get(gene_id)
        if gene is None:
            gene = GeneLevelId.objects.filter(pk=gene_id).first()
            if gene is None:
                continue
            gene_ids[gene_id] = gene
        yield SpliceEventVariant(variant=variant, gene=gene, label=label)


def splice_key(written: str) -> str:
    """ The key a written splice event is matched on - lower-cased, with spaces, hyphens and
        underscores dropped, so 'AR-V7', 'AR V7' and the 'ARV7' an imported c.HGVS arrives as
        (@see ImportedAlleleInfo._tidy_input_value) are one key. Splitting the string into a gene
        and a label instead isn't decidable once the space is gone. """
    return SPLICE_KEY_SEPARATORS.sub("", written or "").lower()


def splice_events_by_key() -> dict[str, SpliceEvent]:
    """ Every written form of a named junction - the canonical 'AR V7' and the 'AR-V7 splice
        variant' a report writes - against the row it names. A naming table of panel junctions, so
        it is small enough to build per call. A label means the same event in every build, so the
        builds' rows collapse onto one key. """

    events: dict[str, SpliceEvent] = {}
    for splice_event in SpliceEvent.objects.all():
        for written in (f"{splice_event.gene_symbol_id} {splice_event.label}", splice_event.display):
            events.setdefault(splice_key(written), splice_event)
    return events


def parse_splice_string(splice_string: str) -> Optional[tuple[str, str, Optional[SpliceEvent]]]:
    """ (gene name, label, the SpliceEvent naming the junction or None) - the gene is returned for a
        resolver to place. A named form is matched whole; the coordinate form names a junction no
        row names, so it is the one form split into its parts. """

    if splice_event := splice_events_by_key().get(splice_key(splice_string)):
        return splice_event.gene_symbol_id, splice_event.label, splice_event
    if m := SPLICE_COORDINATE_PATTERN.match(splice_string or ""):
        return m.group("gene"), _written_coordinate_label(m.group("label")), None
    return None


def _written_coordinate_label(label: str) -> str:
    """ 'chrx_66905968_66914514' as coordinate_label would have written it """
    contig, donor, acceptor = label.split("_")
    return f"{contig[3:] if contig.lower().startswith('chr') else contig}_{donor}_{acceptor}".upper()


def resolve_splice_string(splice_string: str, resolver: GeneLevelNameResolver = None) \
        -> Optional[ResolvedSpliceEvent]:
    """ 'AR V7' -> the identity it is stored under, whose variant_coordinate goes through the VCF
        insert pipeline like any other coordinate.

        A named form resolves against a row we hold. The coordinate form names a junction nothing
        has named, so it resolves only to a Variant that already exists - a lab quoting coordinates
        is pointing at something we loaded, not minting an unnamed junction from free text.
        @see ImportedAlleleInfo for where a classification target comes in this way. """

    parsed = parse_splice_string(splice_string)
    if parsed is None:
        return None

    gene_name, label, splice_event = parsed
    if resolver is None:
        resolver = GeneLevelNameResolver()
    if splice_event is None:
        for splice_event_variant in _find_splice_event_variants(gene_name, label, resolver):
            return ResolvedSpliceEvent(gene=splice_event_variant.gene, label=label)
        return None

    if resolved_gene := resolver.resolve_gene(gene_name, allow_unknown=False):
        return ResolvedSpliceEvent(gene=resolved_gene.gene_level_id, label=label,
                                   splice_event=splice_event)
    return None


def find_splice_events_for_string(splice_string: str, resolver: GeneLevelNameResolver = None) \
        -> list[SpliceEventVariant]:
    """ Lookup only - 'AR V7' finds the event if we have it, and mints nothing if we don't. Search
        runs on whatever a user types, so it must not create identities.

        A splice event has no record of its own, so the Variants are read back through their alt. """

    parsed = parse_splice_string(splice_string)
    if parsed is None:
        return []

    gene_name, label, _splice_event = parsed
    if resolver is None:
        resolver = GeneLevelNameResolver()
    return _find_splice_event_variants(gene_name, label, resolver)


def _find_splice_event_variants(gene_name: str, label: str,
                                resolver: GeneLevelNameResolver) -> list[SpliceEventVariant]:
    """ The splice Variants for a gene name and label, creating nothing - so the GeneLevelIds are
        looked up rather than minted """

    symbols = {resolver.canonical_symbol(n) for n in resolver.split_gene_names(gene_name)}
    genes_by_alt = {}
    for gene in GeneLevelId.objects.filter(symbol_str__in=symbols):
        alt = GeneLevelSymbolicAlt.format(GeneLevelSymbolicAlt.SPLICE, gene.alt_namespace,
                                          gene.pk, label)
        genes_by_alt[alt] = gene
    if not genes_by_alt:
        return []

    variant_qs = Variant.objects.filter(Variant.get_gene_level_q(),
                                        locus__position__in=[g.pk for g in genes_by_alt.values()],
                                        alt__seq__in=list(genes_by_alt)) \
                                .select_related("locus", "alt")
    return [SpliceEventVariant(variant=variant, gene=genes_by_alt[variant.alt.seq], label=label)
            for variant in variant_qs]
