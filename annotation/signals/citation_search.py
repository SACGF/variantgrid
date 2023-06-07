import re

from annotation.models import CitationFetchRequest, Citation
from annotation.models.models_citations import CitationSource, CitationIdNormalized
from snpdb.search import search_receiver, SearchInputInstance, SearchExample


@search_receiver(
    search_type=Citation,
    pattern=re.compile(r"((PMID|PUBMED|PMCID)\s*:\s*(?:PMC)?[0-9]+|NBK[0-9]+)", flags=re.IGNORECASE),
    example=SearchExample(
        note="PMID, PMCID or NBK ID",
        examples=["PMID:25741868", "PMCID:PMC23433", "NBK100238"]
    )
)
def search_citations(search_input: SearchInputInstance):
    try:
        normal_id = CitationIdNormalized.normalize_id(search_input.search_string)

        citation = normal_id.get_or_create()
        CitationFetchRequest.fetch_all_now([normal_id])

        messages = []
        input_string = search_input.search_string
        tidy_input = input_string.replace(' ', '').upper()
        if ':' in tidy_input:
            colon_index = tidy_input.index(':')
            search_prefix = CitationSource.from_legacy_code(tidy_input[:colon_index])
            suffix = tidy_input[colon_index + 1:]
            if citation.source != search_prefix or citation.index != suffix:
                messages.append(f'Normalising "{input_string}" to "{citation.id}"')
        if citation.error:
            messages.append("Could not retrieve citation")

        yield citation, messages
    except ValueError:
        pass
