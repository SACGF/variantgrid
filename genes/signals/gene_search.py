import re

from genes.models import Gene
from genes.models_enums import AnnotationConsortium
from snpdb.search import search_receiver, SearchInputInstance, SearchExample

GENE_PATTERN = re.compile(r"^(?P<prefix>ENSG|GENE\s*(?:ID)?\s*:\s*)(?P<id>\d+)$", re.IGNORECASE)
GENE_VERSION_PATTERN = re.compile(r"^ENSG(?P<id>\d+)\.(?P<version>\d+)$", re.IGNORECASE)

@search_receiver(
    search_type=Gene,
    pattern=GENE_PATTERN,
    example=SearchExample(
        note="Ensemble or RefSeq Gene reference",
        examples=["GeneID:2624", "ENSG00000159216"]
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

#
@search_receiver(
    search_type=Gene,
    pattern=GENE_VERSION_PATTERN,
    example=SearchExample(
        note="Ensemble Gene version",
        examples=["ENSG00000159216.14"]
    )
)
def gene_version_search(search_input: SearchInputInstance):
    """ This is only going to be Ensembl as RefSeq genes don't have versions
        We don't have a gene version page, so we'll just return genes (ignoring version for now)
    """
    gene_id = search_input.match.group('id')
    #  version = search_input.match.group('version')
    gene_id = f"ENSG{gene_id}"
    yield Gene.objects.filter(identifier=gene_id, annotation_consortium=AnnotationConsortium.ENSEMBL)
