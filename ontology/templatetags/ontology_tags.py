from django.template import Library

from ontology.ontology_matching import OntologyMeta

register = Library()


@register.inclusion_tag("ontology/tags/ontology_entry.html")
def ontology(data: OntologyMeta):
    return {"ontology": data}