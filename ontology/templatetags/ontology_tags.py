from typing import Optional

from django.template import Library

from ontology.models import OntologyTerm, OntologyTermRelation, GeneDiseaseClassification
from ontology.ontology_matching import OntologyMatch

register = Library()


@register.inclusion_tag("ontology/tags/ontology_match.html")
def ontology_meta(data: OntologyMatch):
    return {"ontology": data}


@register.inclusion_tag("ontology/tags/ontology_term.html")
def ontology_term(data: OntologyTerm):
    return {"term": data}

@register.inclusion_tag("ontology/tags/ontology_relationship.html")
def ontology_relationship(relationship: OntologyTermRelation, term: OntologyTerm):
    low_quality = False
    quality: Optional[str] = None
    if extra := relationship.extra:
        if strongest := extra.get('strongest_classification'):
            allowed_set = GeneDiseaseClassification.get_above_min(GeneDiseaseClassification.STRONG)
            if strongest not in allowed_set:
                low_quality = True
                quality = strongest

    return {
        "relationship": relationship,
        "low_quality": low_quality,
        "quality": quality,
        "term": term
    }
