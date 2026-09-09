"""
Turning the gene names a fusion caller writes into GeneFusion records.

@see snpdb.gene_level_variants for why a fusion becomes a Variant at all.

A caller writes each side of a fusion as a list of the genes overlapping that breakpoint
(`ROS1;GOPC`, `PPARG/AC016683.6`) using whatever symbol was current when the panel was designed
(`SEPT14` where HGNC now says `SEPTIN14`). Identity needs exactly one gene per side, so a side goes
through alias resolution and then picks its best candidate. Everything the caller wrote is kept
per-observation by the loader, so nothing here is lossy.

Where the caller gave a breakpoint, position decides: it does not depend on which symbol the panel
was designed with, and it is the only evidence that survives a symbol being retired to a different
gene. The name is then the tiebreaker between genes that overlap the same position, and the fallback
where nothing does. Genes found this way are recorded on the FusionGeneId so annotation can reach
them without going back through the symbol (@see annotation.gene_level_annotation).

Names HGNC doesn't carry - clone-based identifiers are routine fusion partners - still get an
identity, via FusionGeneId's local id space, so every call the caller made becomes a Variant.
"""
import re
from dataclasses import dataclass, field
from functools import cached_property
from typing import Optional

from django.db.models import Q

from genes.gene_matching import GeneSymbolMatcher
from genes.gene_overlaps import GeneOverlap, SVGeneOverlapResolver
from genes.models import GeneAnnotationRelease, GeneFusion, FusionGeneId, HGNC, fusion_canonical_str
from genes.models_enums import HGNCStatus
from library.genomics.vcf_enums import GeneLevelSymbolicAlt
from snpdb.clingen_allele import get_variant_allele_for_variant
from snpdb.gene_level_variants import (
    GENE_LEVEL_CONTIG_NAME,
    GENE_LEVEL_REF,
    GENE_LEVEL_SVLEN,
)
from snpdb.models import GenomeBuild, Variant, VariantCoordinate

# Within one cell - a hyphen can't separate, as clone-based identifiers contain them (RP11-458D21.5)
GENE_LIST_SEPARATOR = re.compile(r"[;/]")
# A fusion written as one string, eg by a lab submitting 'BCR::ABL1' as a classification target
FUSION_STRING_SEPARATOR = re.compile(r"::|~|/|--")
# 'chr3:132036420' - what a caller writes as one side's breakpoint
BREAKPOINT = re.compile(r"^\s*(?P<chrom>[^:\s]+)\s*:\s*(?P<position>[0-9,]+)\s*$")


def parse_breakpoint(breakpoint: Optional[str]) -> Optional[tuple[str, int]]:
    """ 'chr3:132036420' -> ('chr3', 132036420). The chromosome is resolved against the build
        later, so 'chr3' and '3' are both fine here """
    if breakpoint and (m := BREAKPOINT.match(breakpoint)):
        return m.group("chrom"), int(m.group("position").replace(",", ""))
    return None


@dataclass
class _GeneCandidate:
    """ The genes at a breakpoint that are all one gene - the same gene in a RefSeq and an Ensembl
        release is two Gene rows, and every release of a build is consulted """
    hgnc_id: Optional[int]
    symbol: Optional[str]
    gene_ids: set[str] = field(default_factory=set)


@dataclass(frozen=True)
class ResolvedFusionGene:
    """ One side of a fusion, resolved to the identity it will be stored under """
    written: str  # The cell exactly as the caller wrote it
    fusion_gene_id: FusionGeneId

    @property
    def resolved_symbol(self) -> str:
        return self.fusion_gene_id.symbol_str

    @property
    def was_renamed(self) -> bool:
        return self.written != self.resolved_symbol


class GeneFusionResolver:
    """ Holds the symbol caches, so build one per file rather than one per row """

    def __init__(self, gene_matcher: GeneSymbolMatcher = None):
        self.gene_matcher = gene_matcher or GeneSymbolMatcher()
        self._overlap_resolvers: dict[int, list[SVGeneOverlapResolver]] = {}

    @staticmethod
    def split_gene_names(cell: str) -> list[str]:
        names = [n.strip() for n in GENE_LIST_SEPARATOR.split(cell or "")]
        return [n for n in names if n]

    @staticmethod
    def split_fusion_string(fusion_string: str) -> Optional[tuple[str, str]]:
        """ 'BCR::ABL1' -> ('BCR', 'ABL1'). A single hyphen is ambiguous with the hyphens inside
            clone-based identifiers, so it only separates when exactly one hyphen is present """
        parts = FUSION_STRING_SEPARATOR.split(fusion_string.strip())
        if len(parts) == 1:
            parts = fusion_string.strip().split("-")
        if len(parts) == 2:
            gene_a, gene_b = (p.strip() for p in parts)
            if gene_a and gene_b:
                return gene_a, gene_b
        return None

    def canonical_symbol(self, name: str) -> str:
        """ The approved symbol where 'name' is an alias (SEPT14 -> SEPTIN14), else the name as given """
        gene_symbol_id, _hgnc = self.resolve_name(name)
        return gene_symbol_id or name

    @cached_property
    def _hgnc_by_previous_symbol(self) -> dict[str, HGNC]:
        """ HGNC's previous symbols, upper-cased - a rename, so the strongest statement that an old
            name and a current one are the same gene """
        return self._hgnc_symbol_lookup("previous_symbols")

    @cached_property
    def _hgnc_by_alias_symbol(self) -> dict[str, HGNC]:
        """ HGNC's alias symbols - a nickname rather than a rename, and routinely shared with an
            unrelated gene, so the last thing consulted """
        return self._hgnc_symbol_lookup("alias_symbols")

    @staticmethod
    def _hgnc_symbol_lookup(field: str) -> dict[str, HGNC]:
        lookup: dict[str, HGNC] = {}
        for hgnc in HGNC.objects.all():
            approved = hgnc.status == HGNCStatus.APPROVED
            for symbol in (getattr(hgnc, field) or "").split(","):
                symbol = symbol.strip().upper()
                if not symbol:
                    continue
                existing = lookup.get(symbol)
                if existing is None or (approved and existing.status != HGNCStatus.APPROVED):
                    lookup[symbol] = hgnc
        return lookup

    @staticmethod
    def _hgnc_for_symbol(gene_symbol_id: str) -> Optional[HGNC]:
        hgnc_qs = HGNC.objects.filter(gene_symbol_id=gene_symbol_id)
        # gene_symbol isn't unique in HGNC (eg MMP21 has multiple entries) so prefer the approved one
        return hgnc_qs.filter(status=HGNCStatus.APPROVED).first() or hgnc_qs.first()

    def resolve_name(self, name: str) -> tuple[Optional[str], Optional[HGNC]]:
        """ (gene symbol, HGNC) - the symbol is the approved one where 'name' was an old name.

            A name that is a current symbol in its own right is taken as written. Otherwise HGNC's
            rename is the strongest evidence two names are one gene, so it outranks an alias: SEPT2
            is SEPTIN2's previous symbol and also, on an unrelated gene, one of SEPTIN6's aliases,
            and GeneSymbolAlias holds a row for each of them. """

        gene_symbol_id, alias_id = self.gene_matcher.get_gene_symbol_id_and_alias_id(name)
        if gene_symbol_id is not None and alias_id is None:
            if hgnc := self._hgnc_for_symbol(gene_symbol_id):
                return gene_symbol_id, hgnc

        uc_name = name.strip().upper()
        if hgnc := self._hgnc_by_previous_symbol.get(uc_name):
            return hgnc.gene_symbol_id, hgnc
        if gene_symbol_id is not None:
            if hgnc := self._hgnc_for_symbol(gene_symbol_id):
                return gene_symbol_id, hgnc

        # A GeneSymbol row for the old name stops the matcher ever reaching its aliases - Ensembl
        # still calls ACP3 'ACPP' - so that hop is taken explicitly here
        if alias_symbol_id := self.gene_matcher.get_alias_gene_symbol_id(name):
            if hgnc := self._hgnc_for_symbol(alias_symbol_id):
                return alias_symbol_id, hgnc
        if hgnc := self._hgnc_by_alias_symbol.get(uc_name):
            return hgnc.gene_symbol_id, hgnc
        return gene_symbol_id, None

    def _gene_overlap_resolvers(self, genome_build: GenomeBuild) -> list[SVGeneOverlapResolver]:
        """ One per GeneAnnotationRelease of the build - a RefSeq release gives the gene its Entrez
            id and an Ensembl one its ENSG, and annotation resolves against whichever it was built
            on. The trees behind these are per contig and built on first use, so a file pays for the
            chromosomes its breakpoints are actually on """
        resolvers = self._overlap_resolvers.get(genome_build.pk)
        if resolvers is None:
            resolvers = [SVGeneOverlapResolver(release)
                         for release in GeneAnnotationRelease.objects.filter(genome_build=genome_build)]
            self._overlap_resolvers[genome_build.pk] = resolvers
        return resolvers

    def _candidates_at_breakpoint(self, breakpoint: Optional[str],
                                  genome_build: Optional[GenomeBuild]) -> list[_GeneCandidate]:
        """ The distinct genes overlapping the breakpoint, grouped by HGNC where they have one """
        position = parse_breakpoint(breakpoint)
        if position is None or genome_build is None:
            return []

        chrom, pos = position
        by_key: dict = {}
        for resolver in self._gene_overlap_resolvers(genome_build):
            for gene_overlap in resolver.get_gene_overlaps(chrom, pos):
                if (key := _candidate_key(gene_overlap)) is None:
                    continue
                candidate = by_key.get(key)
                if candidate is None:
                    candidate = _GeneCandidate(hgnc_id=gene_overlap.hgnc_id, symbol=gene_overlap.symbol)
                    by_key[key] = candidate
                candidate.gene_ids.add(gene_overlap.gene_id)
        return list(by_key.values())

    def _breakpoint_candidate(self, breakpoint: Optional[str], genome_build: Optional[GenomeBuild],
                              resolved_names: list[tuple[Optional[str], Optional[HGNC]]]) \
            -> Optional[_GeneCandidate]:
        """ Which gene at the breakpoint this side is. The name written decides between genes that
            overlap the same position; with nothing to decide on, one candidate stands alone and
            several mean we don't know, so the name is left to answer it """

        candidates = self._candidates_at_breakpoint(breakpoint, genome_build)
        if not candidates:
            return None

        hgnc_ids = {hgnc.pk for _symbol, hgnc in resolved_names if hgnc}
        symbols = {symbol.upper() for symbol, _hgnc in resolved_names if symbol}
        for candidate in candidates:
            if candidate.hgnc_id is not None and candidate.hgnc_id in hgnc_ids:
                return candidate
            if candidate.symbol and candidate.symbol.upper() in symbols:
                return candidate

        if len(candidates) == 1:
            return candidates[0]
        return None

    def _identity_for_candidate(self, candidate: _GeneCandidate) -> FusionGeneId:
        hgnc = HGNC.objects.filter(pk=candidate.hgnc_id).first() if candidate.hgnc_id else None
        if hgnc is not None:
            symbol_str = hgnc.gene_symbol_id
            gene_symbol_id = hgnc.gene_symbol_id
        else:
            symbol_str = candidate.symbol
            gene_symbol_id = candidate.symbol
        fusion_gene_id = FusionGeneId.get_or_create_for_symbol(symbol_str, gene_symbol_id, hgnc)
        fusion_gene_id.genes.add(*candidate.gene_ids)
        return fusion_gene_id

    def resolve_side(self, cell: str, allow_unknown: bool = True, breakpoint: Optional[str] = None,
                     genome_build: Optional[GenomeBuild] = None) -> Optional[ResolvedFusionGene]:
        """ Picks the identity for one side of a fusion.

            The breakpoint decides where the caller gave one and it lands in a gene we know - the
            position does not depend on which symbol the panel was designed with. Failing that an
            HGNC-backed name wins, since that id means the same gene everywhere, then a known gene
            symbol, then the first name as written - a caller naming only clone-based identifiers
            still described a real event.

            allow_unknown=False stops at a known gene symbol, for callers where a name we can't place
            means "this probably isn't a fusion" rather than "this is a fusion of something unusual". """

        names = self.split_gene_names(cell)
        if not names:
            return None

        resolved_names = [self.resolve_name(name) for name in names]
        if candidate := self._breakpoint_candidate(breakpoint, genome_build, resolved_names):
            return ResolvedFusionGene(written=cell, fusion_gene_id=self._identity_for_candidate(candidate))

        first_known_symbol = None
        for gene_symbol_id, hgnc in resolved_names:
            if hgnc is not None:
                fusion_gene_id = FusionGeneId.get_or_create_for_symbol(gene_symbol_id, gene_symbol_id, hgnc)
                return ResolvedFusionGene(written=cell, fusion_gene_id=fusion_gene_id)
            if gene_symbol_id is not None and first_known_symbol is None:
                first_known_symbol = gene_symbol_id

        if first_known_symbol is None and not allow_unknown:
            return None

        symbol_str = first_known_symbol or names[0]
        fusion_gene_id = FusionGeneId.get_or_create_for_symbol(symbol_str, first_known_symbol, None)
        return ResolvedFusionGene(written=cell, fusion_gene_id=fusion_gene_id)

    def resolve_fusion(self, gene_a: Optional[ResolvedFusionGene], gene_b: Optional[ResolvedFusionGene],
                       directionality_known: bool) -> 'ResolvedFusion':
        """ gene_a is the 5' side. Exactly one side may be None - a caller that named neither has
            described nothing to key on. No database writes beyond the FusionGeneIds themselves, so
            the loader can call this before the Variants exist. """

        if gene_a is None and gene_b is None:
            raise ValueError("A fusion needs at least one named gene")

        anchor, partner, is_ordered = _order_partners(gene_a, gene_b, directionality_known)
        return ResolvedFusion(anchor=anchor, partner=partner, is_ordered=is_ordered)


@dataclass(frozen=True)
class ResolvedFusion:
    """ A fusion's identity, before it has a Variant. anchor/partner are the ids that become the
        Locus.position and the alt - @see snpdb.gene_level_variants """
    anchor: FusionGeneId
    partner: Optional[FusionGeneId]
    is_ordered: bool

    @property
    def alt(self) -> str:
        kind = GeneLevelSymbolicAlt.FUSION if self.is_ordered else GeneLevelSymbolicAlt.FUSION_UNORDERED
        if self.partner is not None:
            return GeneLevelSymbolicAlt.format(kind, self.partner.alt_namespace, self.partner.pk)
        return GeneLevelSymbolicAlt.format(kind, None, None)

    @property
    def variant_coordinate(self) -> VariantCoordinate:
        """ What the loader writes as a VCF record, and looks the Variant up by afterwards """
        return VariantCoordinate(chrom=GENE_LEVEL_CONTIG_NAME, position=self.anchor.pk,
                                 ref=GENE_LEVEL_REF, alt=self.alt, svlen=GENE_LEVEL_SVLEN)

    @property
    def canonical_str(self) -> str:
        """ The same string GeneFusion.canonical_str gives, before the Variant exists """
        return fusion_canonical_str(self.anchor, self.partner)


def _candidate_key(gene_overlap: GeneOverlap):
    """ What makes two genes from different releases one candidate. A gene with neither an HGNC nor
        a symbol names nothing an identity could be minted under, so it is no candidate at all """
    if gene_overlap.hgnc_id is not None:
        return "hgnc", gene_overlap.hgnc_id
    if gene_overlap.symbol:
        return "symbol", gene_overlap.symbol.upper()
    return None


def _order_partners(gene_a: Optional[ResolvedFusionGene], gene_b: Optional[ResolvedFusionGene],
                    directionality_known: bool) -> tuple[FusionGeneId, Optional[FusionGeneId], bool]:
    if gene_a is None or gene_b is None:
        # One partner unspecified. Where it's the 5' side that's missing, anchoring the known gene as
        # 5' would assert a direction the caller didn't, so those stay unordered
        known = gene_a or gene_b
        return known.fusion_gene_id, None, gene_b is None and directionality_known

    if directionality_known:
        return gene_a.fusion_gene_id, gene_b.fusion_gene_id, True

    # Unordered anchors on the lower id, so the pair lands on one Variant whichever way it's reported
    anchor, partner = sorted([gene_a.fusion_gene_id, gene_b.fusion_gene_id], key=lambda g: g.pk)
    return anchor, partner, False


def create_gene_fusions_for_variants(variant_qs) -> int:
    """ The GeneFusion rows for gene-level variants an insert pipeline has just created.

        Everything a GeneFusion holds is already in the Variant - the anchor is Locus.position and the
        alt carries the partner and whether a direction was asserted - so this reads the variants
        rather than the file they came from, and one implementation serves every loader.

        :return: how many were created """

    gene_fusions = []
    for variant in variant_qs.filter(Variant.get_gene_level_q(), genefusion__isnull=True) \
                             .select_related("locus", "alt"):
        kind, _namespace, partner_id = GeneLevelSymbolicAlt.parse(variant.alt.seq)
        gene_fusions.append(GeneFusion(variant=variant,
                                       anchor_id=variant.locus.position,
                                       partner_id=partner_id,
                                       is_ordered=kind == GeneLevelSymbolicAlt.FUSION))
    created = GeneFusion.objects.bulk_create(gene_fusions, ignore_conflicts=True)
    return len(created)


def resolve_fusion_string(fusion_string: str, resolver: GeneFusionResolver = None) -> Optional[ResolvedFusion]:
    """ 'BCR::ABL1' -> the identity it will be stored under, whose variant_coordinate goes through the
        VCF insert pipeline like any other coordinate.

        Direction is taken from the order written, which is the convention every fusion nomenclature
        uses. Both sides have to be genes we already know, so an arbitrary hyphenated string doesn't
        mint a fusion. @see ImportedAlleleInfo for where a classification target comes in this way. """

    if resolver is None:
        resolver = GeneFusionResolver()
    if genes := resolver.split_fusion_string(fusion_string):
        gene_a = resolver.resolve_side(genes[0], allow_unknown=False)
        gene_b = resolver.resolve_side(genes[1], allow_unknown=False)
        if gene_a and gene_b:
            return resolver.resolve_fusion(gene_a, gene_b, directionality_known=True)
    return None


def find_gene_fusions_for_string(fusion_string: str, resolver: GeneFusionResolver = None) -> list[GeneFusion]:
    """ Lookup only - 'BCR::ABL1' finds the fusion if we have it, and mints nothing if we don't.

        Search runs on whatever a user types, so it must not create identities. Both orderings come
        back, since someone searching a pair rarely means only one direction, plus the unordered form. """

    if resolver is None:
        resolver = GeneFusionResolver()
    genes = resolver.split_fusion_string(fusion_string)
    if genes is None:
        return []

    gene_ids = []
    for name in genes:
        symbols = {resolver.canonical_symbol(n) for n in resolver.split_gene_names(name)}
        ids = set(FusionGeneId.objects.filter(symbol_str__in=symbols).values_list("pk", flat=True))
        if not ids:
            return []
        gene_ids.append(ids)

    first, second = gene_ids
    q = (Q(anchor__in=first) & Q(partner__in=second)) | (Q(anchor__in=second) & Q(partner__in=first))
    return list(GeneFusion.objects.filter(q).select_related("variant", "anchor", "partner"))


def get_gene_fusion_allele(gene_fusion: GeneFusion, genome_build: GenomeBuild):
    """ Fusion variants sit on a contig every build shares, so one Allele serves them all. ClinGen
        can't register them (no coordinate) - clingen_allele_skip_reason says so - which leaves the
        ordinary 'no ClinGen' path in get_variant_allele_for_variant """
    variant_allele = get_variant_allele_for_variant(genome_build, gene_fusion.variant)
    return variant_allele.allele
