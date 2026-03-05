from seqauto.models import EnrichmentKit
from snpdb.search import search_receiver, SearchInputInstance, SearchExample, HAS_3_ANY


@search_receiver(
    search_type=EnrichmentKit,
    pattern=HAS_3_ANY,
    example=SearchExample(note="3 or more characters of the enrichment kit name")
)
def enrichment_kit_search(search_input: SearchInputInstance):
    yield EnrichmentKit.objects.filter(search_input.q_words())
