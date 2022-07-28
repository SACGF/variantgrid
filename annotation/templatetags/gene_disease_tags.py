from django.template import Library

from ontology.models import OntologySnake, GeneDiseaseClassification, OntologyVersion

register = Library()


@register.inclusion_tag("annotation/tags/gene_disease.html")
def gene_disease(gene_symbol):
    ontology_version = OntologyVersion.latest()
    gene_disease_relations = ontology_version.gene_disease_relations(gene_symbol,
                                                                     min_classification=GeneDiseaseClassification.DISPUTED)
    return {
        "gene_symbol": gene_symbol,
        "gene_disease_relations": gene_disease_relations,
    }
