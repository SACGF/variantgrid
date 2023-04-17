from typing import Any
from django.dispatch import receiver
from genes.models import Gene
from genes.models_enums import AnnotationConsortium
from snpdb.search2 import search_signal, SearchInput, SearchResponse
import re


GENE_PATTERN = re.compile(r"(ENSG|Gene.*:)\d+")


@receiver(search_signal, sender=SearchInput)
def gene_search(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
    """ Symbols have been separated into search_gene_symbol - this returns Gene objects """
    if search_input.matches_pattern(GENE_PATTERN):
        CONSORTIUM_REGEX = {
            r"(ENSG\d+)": AnnotationConsortium.ENSEMBL,
            r"Gene:(\d+)": AnnotationConsortium.REFSEQ,
            r"GeneID:(\d+)": AnnotationConsortium.REFSEQ,
            r"Gene ID:(\d+)": AnnotationConsortium.REFSEQ,
        }

        for c_regex, annotation_consortium in CONSORTIUM_REGEX.items():
            if m := re.match(c_regex, search_input.search_string, re.IGNORECASE):
                gene_id = m.group(1)
                response = SearchResponse(Gene)
                response.extend(Gene.objects.filter(identifier=gene_id, annotation_consortium=annotation_consortium), annotation_consortium=annotation_consortium)
                return response