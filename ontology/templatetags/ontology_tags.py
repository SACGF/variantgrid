from django.template import Library

from ontology.models import OntologyTerm, OntologyRelation
from ontology.ontology_matching import OntologyMatch

register = Library()


@register.inclusion_tag("ontology/tags/ontology_match.html")
def ontology_meta(data: OntologyMatch):
    return {"ontology": data}


@register.inclusion_tag("ontology/tags/ontology_term.html")
def ontology_term(data: OntologyTerm):
    return {"term": data}

@register.inclusion_tag("ontology/tags/ontology_relationship.html")
def ontology_relationship(relationship: OntologyRelation, term: OntologyTerm):
    return {
        "relationship": relationship,
        "term": term
    }