import re

from django.contrib.auth.models import User
from django.dispatch.dispatcher import receiver

from classification.models import ClassificationGroupingSearchTerm, ClassificationGrouping, ClassificationModification
from classification.models.classification_utils import classification_gene_symbol_filter
from classification.signals.classification_search import classification_qs_to_extras
from genes.models import GeneSymbol, GeneSymbolAlias
from library.preview_request import preview_extra_signal, PreviewKeyValue
from snpdb.search import search_receiver, SearchInputInstance, SearchExample, SearchResultMatchStrength

GENE_SYMBOL_PATTERN = re.compile(r"^[a-zA-Z][\da-zA-Z\-\.]+$")


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


@receiver(preview_extra_signal, sender=GeneSymbol)
def gene_symbol_preview_classifications_extra(sender, user: User, obj: GeneSymbol, **kwargs) -> list[PreviewKeyValue]:
    extras = []

    q = classification_gene_symbol_filter(obj)
    cms_qs = ClassificationModification.latest_for_user(user=user).filter(q)
    if cms_qs.exists():
        extras += classification_qs_to_extras(cms_qs)
    return extras
