import operator
from functools import reduce

from django.db.models.query_utils import Q

from annotation.models import PATIENT_ONTOLOGY_TERM_PATH
from genes.models import GeneSymbol, SampleGeneList
from ontology.models import OntologyTerm


def get_sample_ontology_q(ontology_terms_str) -> Q | None:
    sample_patient_ontology_path = f"patient__{PATIENT_ONTOLOGY_TERM_PATH}"
    ontology_filters = []
    for term_name in ontology_terms_str.split(","):
        if ot := OntologyTerm.get_or_stub(term_name):
            ontology_filters.append(Q(**{sample_patient_ontology_path: ot}))
    q = None
    if ontology_filters:
        q = reduce(operator.or_, ontology_filters)
    return q

def get_sample_qc_gene_list_gene_symbol_q(gene_symbol_str) -> Q | None:
    q = None
    if gene_symbol := GeneSymbol.objects.filter(pk=gene_symbol_str).first():
        sample_gene_list_qs = SampleGeneList.objects.filter(gene_list__genelistgenesymbol__gene_symbol=gene_symbol)
        q = Q(samplegenelist__in=sample_gene_list_qs)
    return q
