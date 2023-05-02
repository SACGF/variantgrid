import re
from genes.models import GeneSymbol, GeneSymbolAlias
from snpdb.search import search_receiver, SearchInputInstance, SearchExample, SearchResultMatchStrength

GENE_SYMBOL_PATTERN = re.compile(r"^[a-zA-Z][\da-zA-Z0-9-]+$")


@search_receiver(
    search_type=GeneSymbol,
    pattern=GENE_SYMBOL_PATTERN,
    example=SearchExample(
        note="Full gene symbol identifier",
        examples=["gata2", "BRCA2"]
    ),
    match_strength=SearchResultMatchStrength.ID_MATCH
)
def gene_symbol_alias_search(search_input: SearchInputInstance):
    gene_symbols = GeneSymbol.objects.filter(symbol=search_input.search_string)
    yield gene_symbols

    aliases = GeneSymbolAlias.objects.filter(alias=search_input.search_string).exclude(gene_symbol__in=gene_symbols)
    for alias in aliases:
        yield alias.gene_symbol, f"{alias.alias} is an alias for {alias.gene_symbol}"

