from django.template import Library

from ontology.models import OntologyTerm
from ontology.ontology_matching import OntologyMeta

register = Library()


@register.inclusion_tag("ontology/tags/ontology_meta.html")
def ontology_meta(data: OntologyMeta):
    return {"ontology": data}


@register.inclusion_tag("ontology/tags/ontology_term.html")
def ontology_term(data: OntologyTerm):
    return {"term": data}