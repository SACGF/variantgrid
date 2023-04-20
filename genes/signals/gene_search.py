import re
from genes.models import Gene
from genes.models_enums import AnnotationConsortium
from snpdb.search2 import search_receiver, SearchInputInstance, SearchExample

GENE_PATTERN = re.compile(r"(ENSG|Gene.*:)\d+")


@search_receiver(
    search_type=Gene,
    pattern=GENE_PATTERN,
    example=SearchExample(
        note="Ensemble or RefSeq Gene reference",
        example="GeneID:144"
    )
)
def gene_search(search_input: SearchInputInstance):
    """ Symbols have been separated into search_gene_symbol - this returns Gene objects """
    CONSORTIUM_REGEX = {
        r"(ENSG\d+)": AnnotationConsortium.ENSEMBL,
        r"Gene:(\d+)": AnnotationConsortium.REFSEQ,
        r"GeneID:(\d+)": AnnotationConsortium.REFSEQ,
        r"Gene ID:(\d+)": AnnotationConsortium.REFSEQ,
    }

    for c_regex, annotation_consortium in CONSORTIUM_REGEX.items():
        if m := re.match(c_regex, search_input.search_string, re.IGNORECASE):
            gene_id = m.group(1)
            yield Gene.objects.filter(identifier=gene_id, annotation_consortium=annotation_consortium)