from typing import Any
from django.dispatch import receiver
from genes.models import GeneSymbol, GeneSymbolAlias
from snpdb.search2 import search_signal, SearchInput, SearchResponse
import re


GENE_SYMBOL_PATTERN = re.compile(r"^[a-zA-Z][\da-zA-Z0-9-]+")


@receiver(search_signal, sender=SearchInput)
def gene_symbol_alias_search(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
    if search_input.matches_pattern(GENE_SYMBOL_PATTERN):
        response = SearchResponse(GeneSymbol)
        search_string = search_input.search_string

        gene_symbols = GeneSymbol.objects.filter(symbol=search_string)
        response.extend(gene_symbols)

        aliases = GeneSymbolAlias.objects.filter(alias=search_string).exclude(gene_symbol__in=gene_symbols)
        for alias in aliases:
            response.add(alias, messages=[f"{alias.alias} is an alias for {alias.gene_symbol}"])

        return response
