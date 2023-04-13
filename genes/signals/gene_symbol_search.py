from typing import Type, Any, Optional, List
from django.dispatch import receiver
from genes.models import GeneSymbol, GeneSymbolAlias
from library.preview_request import PreviewData
from snpdb.search2 import SearchResponseRecordAbstract, search_signal, SearchInput, SearchResponse
import re


GENE_SYMBOL_PATTERN = re.compile(r"^[a-zA-Z][\da-zA-Z0-9-]+")


class SearchResponseGeneSymbol(SearchResponseRecordAbstract[GeneSymbol]):

    @classmethod
    def result_class(cls) -> Type:
        return GeneSymbol


@receiver(search_signal, sender=SearchInput)
def gene_symbol_search(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
    # response: SearchResponse[SearchResponseGeneSymbol] = SearchResponse(SearchResponseGeneSymbol)
    # response.mark_valid_search()
    # response.extend(GeneSymbol.objects.filter(symbol="BRCA2"))
    # return response

    if search_input.matches_pattern(GENE_SYMBOL_PATTERN):
        response: SearchResponse[SearchResponseGeneSymbol] = SearchResponse(SearchResponseGeneSymbol)
        response.mark_valid_search()
        response.extend(GeneSymbol.objects.filter(symbol=search_input.search_string.upper()))
        return response


class SearchResponseGeneSymbolAlias(SearchResponseRecordAbstract[GeneSymbolAlias]):

    @classmethod
    def result_class(cls) -> Type:
        return GeneSymbolAlias

    @classmethod
    def category(cls) -> str:
        return SearchResponseGeneSymbol.category()

    @property
    def preview(self) -> PreviewData:
        return self.record.gene_symbol.preview

    @property
    def messages(self) -> Optional[List[str]]:
        return [f"{self.record.alias} is an alias for {self.record.gene_symbol}"]


@receiver(search_signal, sender=SearchInput)
def gene_symbol_alias_search(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
    if search_input.matches_pattern(GENE_SYMBOL_PATTERN):
        response: SearchResponse[SearchResponseGeneSymbol] = SearchResponse(SearchResponseGeneSymbolAlias)
        search_string = search_input.search_string

        gene_symbols = GeneSymbol.objects.filter(symbol=search_string)
        aliases = GeneSymbolAlias.objects.filter(alias=search_string).exclude(gene_symbol__in=gene_symbols)

        response.extend(aliases)
        return response
