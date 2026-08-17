from seqauto.models import EnrichmentKit
from snpdb.search import HAS_3_ANY, SearchExample, SearchInputInstance, search_receiver


@search_receiver(
    search_type=EnrichmentKit,
    pattern=HAS_3_ANY,
    example=SearchExample(note="3 or more characters of the enrichment kit name")
)
def enrichment_kit_search(search_input: SearchInputInstance):
    yield EnrichmentKit.objects.filter(search_input.q_words())
