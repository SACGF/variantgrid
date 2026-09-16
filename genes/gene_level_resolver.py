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
"""
import re
from dataclasses import dataclass
from functools import cached_property
from typing import Optional

from genes.gene_matching import GeneSymbolMatcher
from genes.models import HGNC, GeneLevelId
from genes.models_enums import HGNCStatus

# Within one cell - a hyphen can't separate, as clone-based identifiers contain them (RP11-458D21.5)
GENE_LIST_SEPARATOR = re.compile(r"[;/]")


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
