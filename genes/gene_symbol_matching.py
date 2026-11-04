"""
The Gene Symbol Matching page (#2637): the user-facing statement of GeneSymbol.resolve_names. For each name it
says which rule decided the result - taken as written, a rename, a single-gene alias, an ambiguous alias, or
nothing known - over the worked examples the page shows and the name a user typed into it.

Entry point: explain_gene_names.
"""
from collections.abc import Iterable
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

from genes.models.models_gene import HGNC, GeneNameResolution, GeneSymbol, GeneSymbolAlias

# Names the page resolves live, and what each one illustrates
GENE_NAME_EXAMPLES: list[tuple[str, str]] = [
    ("BRCA2", "A current symbol is taken as written."),
    ("SEPT2", "SEPTIN2's previous symbol, and also one of SEPTIN6's aliases - the rename wins."),
    ("ACPP", "An old symbol that still has its own row here (Ensembl kept writing it) moves to its rename, ACP3."),
    ("AIP", "An approved symbol that another gene (AURKAIP1) lists as an alias stays itself."),
    ("A20", "An alias of several genes is refused rather than one of them being picked."),
]


class GeneNameRule(Enum):
    """ Why a name resolved the way it did - the branches of genes.models.models_gene.resolve_gene_name """
    CURRENT_SYMBOL = "A current HGNC symbol, so it is taken as written"
    PREVIOUS_SYMBOL = "A previous symbol - HGNC renamed the gene, so the current symbol is used"
    ALIAS = "An alias of a single gene, so that gene is used"
    AMBIGUOUS = "An alias of several genes, so it is not matched - use the gene's current symbol"
    SYMBOL_WITHOUT_HGNC = "A symbol in this database with no HGNC record, kept as written"
    ALIAS_WITHOUT_HGNC = "An alias of a symbol HGNC no longer carries, kept as that symbol"
    UNKNOWN = "Not a symbol or alias this database knows"


@dataclass(frozen=True)
class GeneNameExplanation:
    name: str
    rule: GeneNameRule
    gene_symbol_id: Optional[str] = None
    hgnc: Optional[HGNC] = None
    alias: Optional[GeneSymbolAlias] = None  # The GeneSymbolAlias row the name was followed through
    candidates: list[str] = field(default_factory=list)  # The genes an ambiguous alias names

    @property
    def matched(self) -> bool:
        return self.gene_symbol_id is not None


def _lists_symbol(hgnc: HGNC, field_name: str, name: str) -> bool:
    symbols = (getattr(hgnc, field_name) or "").split(",")
    return name.upper() in {symbol.strip().upper() for symbol in symbols}


def _rule_for(name: str, resolution: GeneNameResolution, alias: Optional[GeneSymbolAlias]) -> GeneNameRule:
    if resolution.gene_symbol_id is None:
        return GeneNameRule.UNKNOWN
    if alias is None:
        return GeneNameRule.CURRENT_SYMBOL if resolution.hgnc else GeneNameRule.SYMBOL_WITHOUT_HGNC
    if resolution.hgnc is None:
        return GeneNameRule.ALIAS_WITHOUT_HGNC
    if _lists_symbol(resolution.hgnc, "previous_symbols", name):
        return GeneNameRule.PREVIOUS_SYMBOL
    return GeneNameRule.ALIAS


def explain_gene_names(names: Iterable[str]) -> dict[str, GeneNameExplanation]:
    """ Each name as passed -> how GeneSymbol.resolve_names treats it and why. One resolve plus one alias query """
    resolutions = GeneSymbol.resolve_names(names)
    alias_ids = {resolution.gene_symbol_alias_id for resolution, _ambiguous in resolutions.values()
                 if resolution.gene_symbol_alias_id}
    aliases = GeneSymbolAlias.objects.in_bulk(alias_ids)

    explanations = {}
    for name, (resolution, ambiguous) in resolutions.items():
        if ambiguous:
            explanations[name] = GeneNameExplanation(name=name, rule=GeneNameRule.AMBIGUOUS,
                                                     candidates=ambiguous.gene_symbol_ids)
            continue
        alias = aliases.get(resolution.gene_symbol_alias_id)
        explanations[name] = GeneNameExplanation(name=name, rule=_rule_for(name, resolution, alias),
                                                 gene_symbol_id=resolution.gene_symbol_id,
                                                 hgnc=resolution.hgnc, alias=alias)
    return explanations
