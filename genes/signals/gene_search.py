import re

from genes.models import Gene
from genes.models_enums import AnnotationConsortium
from snpdb.search import search_receiver, SearchInputInstance, SearchExample

GENE_PATTERN = re.compile(r"^(?P<prefix>ENSG|GENE\s*(?:ID)?\s*:\s*)(?P<id>\d+)$", re.IGNORECASE)


@search_receiver(
    search_type=Gene,
    pattern=GENE_PATTERN,
    example=SearchExample(
        note="Ensemble or RefSeq Gene reference",
        examples=["GeneID:2624"]
    )
)
def gene_search(search_input: SearchInputInstance):
    prefix = search_input.match.group('prefix')
    gene_id = search_input.match.group('id')
    consortium = AnnotationConsortium.REFSEQ
    if prefix.upper().startswith("ENSG"):
        consortium = AnnotationConsortium.ENSEMBL
        gene_id = f"ENSG{gene_id}"

    yield Gene.objects.filter(identifier=gene_id, annotation_consortium=consortium)