"""
Turning a gene name a caller wrote into the GeneLevelId it is stored under.

Shared by both kinds of gene-level event: a fusion resolves one of these per side and can let the
breakpoint decide (@see genes.gene_fusions), a whole-gene copy number call resolves the one gene its
segment names (@see genes.gene_copy_number). @see snpdb.gene_level_variants for why either is a
Variant at all.

A caller writes a gene as a list of the symbols at that place (`ROS1;GOPC`, `PPARG/AC016683.6`)
using whatever symbol was current when the panel was designed (`SEPT14` where HGNC now says
`SEPTIN14`, `MYCL1` where it says `MYCL`). Identity needs exactly one gene, so a cell goes through
alias resolution and then picks its best candidate.

Names HGNC doesn't carry - clone-based identifiers are routine fusion partners - still get an
identity, via GeneLevelId's local id space, so every call the caller made becomes a Variant.

Where a caller gives a position, GenePositionResolver finds the gene there instead - a fusion
breakpoint, or a splice junction whose caller names no gene at all (@see genes.gene_splice).
"""
import re
from dataclasses import dataclass, field
from functools import cached_property
from typing import Any, Optional

from genes.gene_matching import GeneSymbolMatcher
from genes.gene_overlaps import GeneOverlap, SVGeneOverlapResolver
from genes.models import HGNC, GeneAnnotationRelease, GeneLevelId
from genes.models_enums import HGNCStatus
from snpdb.models import GenomeBuild

# Within one cell - a hyphen can't separate, as clone-based identifiers contain them (RP11-458D21.5)
GENE_LIST_SEPARATOR = re.compile(r"[;/]")


@dataclass(frozen=True)
class GeneLevelResolution:
    """ What a resolve_*_string call answers: the identity the string names, or the reason the kind
        of event it was recognised as refused it.

        Neither set means the string is not that kind of event at all, so the next resolver gets a
        turn; a reason means it was that kind and did not validate, which is what a record's message
        says (@see classification.models.ImportedAlleleInfo). """
    resolved: Optional[Any] = None
    reason: Optional[str] = None

    def __bool__(self) -> bool:
        return self.resolved is not None

    @property
    def recognised(self) -> bool:
        """ The string is this kind of gene-level event, resolved or refused """
        return self.resolved is not None or self.reason is not None

    @staticmethod
    def not_applicable() -> 'GeneLevelResolution':
        return GeneLevelResolution()

    @staticmethod
    def identity(resolved: Any) -> 'GeneLevelResolution':
        return GeneLevelResolution(resolved=resolved)

    @staticmethod
    def refused(reason: str) -> 'GeneLevelResolution':
        return GeneLevelResolution(reason=reason)


def unknown_gene_reason(name: str) -> str:
    """ The reason every gene-level resolver gives for a name that is no symbol we hold - a lab
        typo'd partner, so the message names the half that failed """
    return f"gene '{name}' is not a symbol we know"


@dataclass(frozen=True)
class ResolvedGeneLevelGene:
    """ One gene of a gene-level event, resolved to the identity it will be stored under """
    written: str  # The cell exactly as the caller wrote it
    gene_level_id: GeneLevelId

    @property
    def resolved_symbol(self) -> str:
        return self.gene_level_id.symbol_str

    @property
    def was_renamed(self) -> bool:
        return self.written != self.resolved_symbol


class GeneLevelNameResolver:
    """ Holds the symbol caches, so build one per file rather than one per row """

    def __init__(self, gene_matcher: GeneSymbolMatcher = None):
        self.gene_matcher = gene_matcher or GeneSymbolMatcher()

    @staticmethod
    def split_gene_names(cell: str) -> list[str]:
        names = [n.strip() for n in GENE_LIST_SEPARATOR.split(cell or "")]
        return [n for n in names if n]

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

    def resolve_gene(self, cell: str, allow_unknown: bool = True) -> Optional[ResolvedGeneLevelGene]:
        """ Picks the identity for one gene, on the name alone.

            An HGNC-backed name wins, since that id means the same gene everywhere, then a known gene
            symbol, then the first name as written - a caller naming only clone-based identifiers
            still described a real event.

            allow_unknown=False stops at a known gene symbol, for callers where a name we can't place
            means "this probably isn't what it looks like" rather than "this is an unusual gene". """

        names = self.split_gene_names(cell)
        if not names:
            return None
        return self._resolve_from_names(cell, names, [self.resolve_name(name) for name in names],
                                        allow_unknown)

    def _resolve_from_names(self, cell: str, names: list[str],
                            resolved_names: list[tuple[Optional[str], Optional[HGNC]]],
                            allow_unknown: bool) -> Optional[ResolvedGeneLevelGene]:
        first_known_symbol = None
        for gene_symbol_id, hgnc in resolved_names:
            if hgnc is not None:
                gene_level_id = GeneLevelId.get_or_create_for_symbol(gene_symbol_id, gene_symbol_id, hgnc)
                return ResolvedGeneLevelGene(written=cell, gene_level_id=gene_level_id)
            if gene_symbol_id is not None and first_known_symbol is None:
                first_known_symbol = gene_symbol_id

        if first_known_symbol is None and not allow_unknown:
            return None

        symbol_str = first_known_symbol or names[0]
        gene_level_id = GeneLevelId.get_or_create_for_symbol(symbol_str, first_known_symbol, None)
        return ResolvedGeneLevelGene(written=cell, gene_level_id=gene_level_id)


@dataclass
class GeneCandidate:
    """ The genes at a position that are all one gene - the same gene in a RefSeq and an Ensembl
        release is two Gene rows, and every release of a build is consulted """
    hgnc_id: Optional[int]
    symbol: Optional[str]
    gene_ids: set[str] = field(default_factory=set)


def _candidate_key(gene_overlap: GeneOverlap):
    """ What makes two genes from different releases one candidate. A gene with neither an HGNC nor
        a symbol names nothing an identity could be minted under, so it is no candidate at all """
    if gene_overlap.hgnc_id is not None:
        return "hgnc", gene_overlap.hgnc_id
    if gene_overlap.symbol:
        return "symbol", gene_overlap.symbol.upper()
    return None


class GenePositionResolver(GeneLevelNameResolver):
    """ The name resolver plus the genes at a position. Holds the per-contig transcript trees as well
        as the symbol caches, so build one per file rather than one per row """

    def __init__(self, gene_matcher: GeneSymbolMatcher = None):
        super().__init__(gene_matcher)
        self._overlap_resolvers: dict[int, list[SVGeneOverlapResolver]] = {}

    def _gene_overlap_resolvers(self, genome_build: GenomeBuild) -> list[SVGeneOverlapResolver]:
        """ One per GeneAnnotationRelease of the build - a RefSeq release gives the gene its Entrez
            id and an Ensembl one its ENSG, and annotation resolves against whichever it was built
            on. The trees behind these are per contig and built on first use, so a file pays for the
            chromosomes its positions are actually on """
        resolvers = self._overlap_resolvers.get(genome_build.pk)
        if resolvers is None:
            resolvers = [SVGeneOverlapResolver(release)
                         for release in GeneAnnotationRelease.objects.filter(genome_build=genome_build)]
            self._overlap_resolvers[genome_build.pk] = resolvers
        return resolvers

    def candidates_at(self, genome_build: Optional[GenomeBuild], chrom: str, position: int) -> list[GeneCandidate]:
        """ The distinct genes overlapping the position, grouped by HGNC where they have one """
        if genome_build is None:
            return []

        by_key: dict = {}
        for resolver in self._gene_overlap_resolvers(genome_build):
            for gene_overlap in resolver.get_gene_overlaps(chrom, position):
                if (key := _candidate_key(gene_overlap)) is None:
                    continue
                candidate = by_key.get(key)
                if candidate is None:
                    candidate = GeneCandidate(hgnc_id=gene_overlap.hgnc_id, symbol=gene_overlap.symbol)
                    by_key[key] = candidate
                candidate.gene_ids.add(gene_overlap.gene_id)
        return list(by_key.values())

    def identity_for_candidate(self, candidate: GeneCandidate) -> GeneLevelId:
        """ Genes found by position are recorded on the GeneLevelId so annotation can reach them
            without going back through the symbol (@see annotation.gene_level_annotation) """
        hgnc = HGNC.objects.filter(pk=candidate.hgnc_id).first() if candidate.hgnc_id else None
        if hgnc is not None:
            symbol_str = hgnc.gene_symbol_id
            gene_symbol_id = hgnc.gene_symbol_id
        else:
            symbol_str = candidate.symbol
            gene_symbol_id = candidate.symbol
        gene_level_id = GeneLevelId.get_or_create_for_symbol(symbol_str, gene_symbol_id, hgnc)
        gene_level_id.genes.add(*candidate.gene_ids)
        return gene_level_id
