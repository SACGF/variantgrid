"""
Splice calls - "AR-V7", "MET exon 14 skipping" - as gene-level variants.

@see snpdb.gene_level_variants for why one of these is a Variant, and genes.gene_level_resolver for
turning the caller's gene name into the identity it is stored under.

Identity is the gene plus the junction's label, so two events in one gene are two variants. The
label is the lab's own name for the junction, canonicalised (canonical_splice_label): lower-case
tokens joined by underscores - v_7, v_iii, exon_14_skipping, or, for a junction named by its
breakpoints, grch37_x_66905968_66914514. Every written form of one junction canonicalises to one
label, so duplicates are prevented by canonicalisation rather than by a table, and a name nobody
registered ahead of time still mints its Variant.

display_splice_label formats a label back for a human (AR-V7 splice, EGFRvIVa splice, MET exon 14 skipping).
genes.models.models_splice_event.SpliceEvent is consulted only where a row's own wording should win
- the junctions the TSO 500 panel reports - and is never asked whether a name is real.

The alt is a Sequence, so the label on it is upper-cased (<SPLICE:HGNC:7029:EXON_14_SKIPPING>) -
that is the storage form only, and GeneLevelSymbolicAlt.parse lowers it back to the canonical label.

A splice event has no record of its own the way a fusion has a GeneFusion - the alt carries the gene
and the label, and the caller's breakpoints ride along in the VCF's INFO - so the objects here are
resolved identities and a read-only view of a Variant, not models.
"""
import re
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from typing import Optional

from genes.gene_level_resolver import (
    GeneLevelNameResolver,
    GeneLevelResolution,
    unknown_gene_reason,
)
from genes.models import GeneLevelId
from genes.models.models_splice_event import SpliceEvent
from library.genomics.vcf_enums import GeneLevelSymbolicAlt
from snpdb.gene_level_variants import (
    GENE_LEVEL_CONTIG_NAME,
    GENE_LEVEL_REF,
    GENE_LEVEL_SVLEN,
)
from snpdb.models import Contig, GenomeBuild, Variant, VariantCoordinate

# The contig of a junction named by its breakpoints, spelled out rather than left as \w+ so
# 'ARX_66905968_66914514' splits after the gene
CONTIG_TOKEN = r"(?:chr)?(?:[0-9]{1,2}|MT|M|X|Y)"
# Inside a written label - 'AR-V7', 'MET ex14 skipping', and the canonical 'v_7' itself
_GAP = r"[-_\s]*"
# The label shapes we accept, as a lab writes them. Each shape's parts are captured, since
# canonicalising is reading them back out in one order
_NUMBERED_LABEL = rf"v{_GAP}(?P<number>[0-9]+)"
_ROMAN_LABEL = rf"v{_GAP}(?P<roman>[IVX]+)(?P<roman_suffix>[a-z]?)"
_EXON_SKIPPING_LABEL = rf"ex(?:on)?{_GAP}(?P<exon>[0-9]+){_GAP}skip(?:ping)?"
_COORDINATE_LABEL = rf"(?P<contig>{CONTIG_TOKEN})[_:](?P<donor>[0-9]+)[_-](?P<acceptor>[0-9]+)"
SPLICE_LABEL_PATTERN = rf"(?:{_NUMBERED_LABEL}|{_ROMAN_LABEL}|{_EXON_SKIPPING_LABEL}|{_COORDINATE_LABEL})"
# 'splice variant' - what a report adds after the label, and a lab leaves on
_SPLICE_SUFFIX = r"(?:\s*splice(?:\s*variant)?)?"
# 'AR V7', 'AR-V7 splice variant', 'MET exon 14 skipping'. The space between gene and label is
# optional ('ARV7', and the classifications imported before gene-level values kept their spaces);
# a gene that has to resolve is what keeps this from claiming ordinary words.
SPLICE_STRING_PATTERN = re.compile(
    rf"^\s*(?P<gene>[A-Za-z0-9.]+?)\s*[-_]?\s*{SPLICE_LABEL_PATTERN}{_SPLICE_SUFFIX}\s*$",
    re.IGNORECASE)
# The same label on its own, for a label read off an alt or handed in by a test
SPLICE_LABEL_ONLY_PATTERN = re.compile(rf"^\s*[-_]?\s*{SPLICE_LABEL_PATTERN}{_SPLICE_SUFFIX}\s*$",
                                       re.IGNORECASE)
# A label already in canonical form - recognised as it stands, so canonicalising is idempotent and a
# coordinate label (whose build we could not work out again) round trips
CANONICAL_COORDINATE_PATTERN = re.compile(
    r"^(?P<build>.+)_(?P<contig>[0-9]{1,2}|mt|m|x|y)_(?P<donor>[0-9]+)_(?P<acceptor>[0-9]+)$")
CANONICAL_NUMBERED_PATTERN = re.compile(r"^v_(?P<number>[0-9]+)$")
CANONICAL_ROMAN_PATTERN = re.compile(r"^v_(?P<roman>[ivx]+)(?P<roman_suffix>[a-z]?)$")
CANONICAL_EXON_SKIPPING_PATTERN = re.compile(r"^exon_(?P<exon>[0-9]+)_skipping$")
_CANONICAL_PATTERNS = (CANONICAL_NUMBERED_PATTERN, CANONICAL_ROMAN_PATTERN,
                       CANONICAL_EXON_SKIPPING_PATTERN, CANONICAL_COORDINATE_PATTERN)
# What display_splice_label ends a string with, so a label alone says what it is wherever it turns
# up - an annotation column, a log line. Exon skipping already says so and takes none
SPLICE_DISPLAY_SUFFIX = " splice"


def genome_build_token(genome_build: GenomeBuild) -> str:
    """ The build as a label token - 'grch37', 't2t_chm13v2_0'. A junction named by its breakpoints
        means different junctions in different builds, so the build is part of its identity """
    return re.sub(r"[^A-Za-z0-9]+", "_", genome_build.name).lower()


def canonical_splice_label(written_label: str, genome_build: GenomeBuild = None) -> Optional[str]:
    """ 'V7', '-V7 splice variant' -> 'v_7'; 'Exon14Skipping' -> 'exon_14_skipping'. The one place a
        junction's name becomes the label it is stored under, so every written form of one junction
        is one Variant. A label already canonical comes back as it is.

        The breakpoint form takes the build it was written in, which is part of what it names - so
        it needs one, and gives None without it. """

    if not written_label:
        return None
    lowered = written_label.strip().lower()
    if any(pattern.match(lowered) for pattern in _CANONICAL_PATTERNS):
        return lowered
    if m := SPLICE_LABEL_ONLY_PATTERN.match(written_label):
        return _canonical_from_match(m, genome_build)
    return None


def _canonical_from_match(m: re.Match, genome_build: Optional[GenomeBuild]) -> Optional[str]:
    """ The canonical label for a matched written form - the tokens, lower-cased and underscored """
    if number := m.group("number"):
        return f"v_{int(number)}"
    if roman := m.group("roman"):
        return f"v_{roman.lower()}{m.group('roman_suffix').lower()}"
    if exon := m.group("exon"):
        return f"exon_{int(exon)}_skipping"
    if contig := m.group("contig"):
        if genome_build is None:
            return None
        contig = re.sub(r"^chr", "", contig, flags=re.IGNORECASE).lower()
        return f"{genome_build_token(genome_build)}_{contig}_{int(m.group('donor'))}_{int(m.group('acceptor'))}"
    return None


def display_splice_label(gene: GeneLevelId, label: str) -> str:
    """ 'v_iva' -> 'EGFRvIVa splice' - the canonical label written the way the literature writes it
        plus SPLICE_DISPLAY_SUFFIX, which is the form to display and to send anywhere off this
        deployment (@see GeneLevelId). A SpliceEvent row's own display wins where we have one
        (@see SpliceEventVariant.display). """

    symbol = gene.symbol_str
    if m := CANONICAL_NUMBERED_PATTERN.match(label or ""):
        return f"{symbol}-V{m.group('number')}{SPLICE_DISPLAY_SUFFIX}"
    if m := CANONICAL_ROMAN_PATTERN.match(label or ""):
        return f"{symbol}v{m.group('roman').upper()}{m.group('roman_suffix')}{SPLICE_DISPLAY_SUFFIX}"
    if m := CANONICAL_EXON_SKIPPING_PATTERN.match(label or ""):
        return f"{symbol} exon {m.group('exon')} skipping"
    if m := CANONICAL_COORDINATE_PATTERN.match(label or ""):
        build = _genome_build_display(m.group("build"))
        return f"{symbol} {build} {m.group('contig').upper()}:{m.group('donor')}-{m.group('acceptor')}{SPLICE_DISPLAY_SUFFIX}"
    return f"{symbol} {label}{SPLICE_DISPLAY_SUFFIX}"


def _genome_build_display(build_token: str) -> str:
    """ The build's own name for a label token. GenomeBuild's manager caches the table, so this is
        not a query per label """
    for genome_build in GenomeBuild.objects.all():
        if genome_build_token(genome_build) == build_token:
            return genome_build.name
    return build_token


def coordinate_label(genome_build: GenomeBuild, contig: Contig, donor: int, acceptor: int) -> str:
    """ The label for a junction the lab has no name for - its breakpoints in the build they were
        called in, which reads as raw coordinates wherever the label is printed """
    contig_name = contig.name.replace("chr", "")
    return f"{genome_build_token(genome_build)}_{contig_name.lower()}_{donor}_{acceptor}"


@dataclass(frozen=True)
class ResolvedSpliceEvent:
    """ A splice event's identity, before it has a Variant. The gene is the Locus.position and the
        alt carries it again plus the canonical label - @see snpdb.gene_level_variants """
    gene: GeneLevelId
    label: str
    splice_event: Optional[SpliceEvent] = None  # The junction's own wording, where a row has one

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
        return display_splice_label(self.gene, self.label)

    @property
    def display(self) -> str:
        """ What a report writes - the row's own wording where we have one, else the canonical form """
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
        """ The identity a caller's row is stored under. A SpliceEvent for the junction's
            coordinates gives it the label a classification for the same junction arrives under;
            otherwise the breakpoints themselves are the label """

        resolved_gene = self.name_resolver.resolve_gene(gene_name)
        if resolved_gene is None:
            return None

        contig = self.genome_build.chrom_contig_mappings.get(chrom)
        if contig is None:
            return None

        splice_event = self.get_splice_event(contig, donor, acceptor)
        if splice_event:
            label = splice_event.label
        else:
            label = coordinate_label(self.genome_build, contig, donor, acceptor)
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
    def variant_id(self) -> int:
        """ Named as the FK a GeneFusion or GeneCopyNumberEvent has, so anything walking gene-level
            events reads all three the same way """
        return self.variant.pk

    @property
    def gene_level_ids(self) -> list[GeneLevelId]:
        """ The genes this event is about, so gene lists and annotation read it the way they read a
            fusion's partners """
        return [self.gene]

    @property
    def canonical_str(self) -> str:
        return display_splice_label(self.gene, self.label)

    @property
    def splice_event(self) -> Optional[SpliceEvent]:
        """ The junction's own wording, where the TSO 500 panel reports it. Keyed on (gene symbol,
            label) rather than the coordinates, which the Variant does not carry - the pair is
            unique per build, and a label means the same event in every build """
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


def parse_splice_string(splice_string: str, genome_build: GenomeBuild = None) -> Optional[tuple[str, str]]:
    """ 'EGFRvIVa' -> ('EGFR', 'v_iva') - the gene name as written, for a resolver to place, and the
        canonical label. No table is consulted: the lab's label is the identity """

    if m := SPLICE_STRING_PATTERN.match(splice_string or ""):
        if label := _canonical_from_match(m, genome_build):
            return m.group("gene"), label
    return None


def resolve_splice_string(splice_string: str, genome_build: GenomeBuild = None,
                          resolver: GeneLevelNameResolver = None) -> GeneLevelResolution:
    """ 'AR V7' -> the identity it is stored under, whose variant_coordinate goes through the VCF
        insert pipeline like any other coordinate.

        A splice string goes through the stages an HGVS does - canonicalise, validate, coordinate,
        match - and mints its Variant whenever it validates, named or not. The gene has to be one we
        already know, and a junction named by its breakpoints has to name a contig of the build it
        was imported in. @see ImportedAlleleInfo for where a classification target comes in this way. """

    if not SPLICE_STRING_PATTERN.match(splice_string or ""):
        return GeneLevelResolution.not_applicable()

    parsed = parse_splice_string(splice_string, genome_build)
    if parsed is None:
        return GeneLevelResolution.refused(
            f"'{splice_string}' names a junction by its breakpoints, which needs the build it was called in")

    gene_name, label = parsed
    if resolver is None:
        resolver = GeneLevelNameResolver()
    resolved_gene = resolver.resolve_gene(gene_name, allow_unknown=False)
    if resolved_gene is None:
        return GeneLevelResolution.refused(unknown_gene_reason(gene_name))

    if genome_build and (m := CANONICAL_COORDINATE_PATTERN.match(label)):
        contig_name = m.group("contig").upper()
        if genome_build.chrom_contig_mappings.get(contig_name) is None:
            return GeneLevelResolution.refused(f"'{contig_name}' is not a contig of {genome_build}")

    return GeneLevelResolution.identity(ResolvedSpliceEvent(gene=resolved_gene.gene_level_id, label=label))


def find_splice_events_for_string(splice_string: str, genome_build: GenomeBuild = None,
                                  resolver: GeneLevelNameResolver = None) -> list[SpliceEventVariant]:
    """ Lookup only - 'AR V7' finds the event if we have it, and mints nothing if we don't. Search
        runs on whatever a user types, so it must not create identities.

        A splice event has no record of its own, so the Variants are read back through their alt. """

    parsed = parse_splice_string(splice_string, genome_build)
    if parsed is None:
        return []

    gene_name, label = parsed
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
