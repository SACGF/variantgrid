import uuid

from django.db.models.aggregates import Count
from django.template import Library

from genes.models import GeneSymbol
from ontology.models import OntologyTerm, OntologyService

GENES_COLOR = "#99CD83"
HPO_COLOR = "#FF8C00"
OMIM_COLOR = "#22529e"
MATCHES_GRAPH_TEMPLATE = "patients/tags/patient_matches_graph.html"

register = Library()


def match_graph(title, qs, field_label, color, graph_width, graph_height, max_records, click_handler=None, patient_ids=None):
    patient_path = "textphenotypematch__text_phenotype__textphenotypesentence__phenotype_description__patienttextphenotype__patient"
    filter_kwargs = {patient_path + "__isnull": False}
    if patient_ids:
        filter_kwargs[patient_path + "__in"] = patient_ids

    qs = qs.filter(**filter_kwargs).annotate(count=Count(patient_path, distinct=True))
    qs = qs.order_by("-count").values_list(field_label, "count")
    if max_records:
        qs = qs[:max_records]

    labels = []
    counts = []
    for n, count in qs:
        labels.append(n)
        counts.append(count)

    return {'title': "%s<br>Patient Count" % title,
            'uuid': uuid.uuid4(),
            'graph_width': graph_width,
            'graph_height': graph_height,
            'color': color,
            'labels': labels,
            'counts': counts,
            'click_handler': click_handler}


@register.inclusion_tag(MATCHES_GRAPH_TEMPLATE)
def patient_phenotypes_graph(graph_width=512, graph_height=384, max_records=20, click_handler=None, patient_ids=None):
    qs = OntologyTerm.objects.filter(ontology_service=OntologyService.HPO)
    return match_graph("Human Phenotype Ontology", qs, "name", HPO_COLOR, graph_width=graph_width, graph_height=graph_height, max_records=max_records, click_handler=click_handler, patient_ids=patient_ids)


@register.inclusion_tag(MATCHES_GRAPH_TEMPLATE)
def patient_omim_graph(graph_width=512, graph_height=384, max_records=20, click_handler=None, patient_ids=None):
    qs = OntologyTerm.objects.filter(ontology_service=OntologyService.OMIM)
    return match_graph("OMIM", qs, "name", OMIM_COLOR, graph_width=graph_width, graph_height=graph_height, max_records=max_records, click_handler=click_handler, patient_ids=patient_ids)


@register.inclusion_tag(MATCHES_GRAPH_TEMPLATE)
def patient_genes_graph(graph_width=512, graph_height=384, max_records=20, click_handler=None, patient_ids=None):
    qs = GeneSymbol.objects.all()
    return match_graph("Genes", qs, "symbol", GENES_COLOR, graph_width=graph_width, graph_height=graph_height, max_records=max_records, click_handler=click_handler, patient_ids=patient_ids)
