from collections import Counter, defaultdict
from django.template import Library

from ontology.models import GeneDiseaseClassification, OntologyVersion

register = Library()


@register.inclusion_tag("annotation/tags/gene_disease.html")
def gene_disease(gene_symbol):
    ontology_version = OntologyVersion.latest()
    try:
        gene_disease_relations = ontology_version.gene_disease_relations(gene_symbol,
                                                                         min_classification=GeneDiseaseClassification.DISPUTED)
    except ValueError:  # No HGNC for symbol
        gene_disease_relations = None

    return {
        "gene_symbol": gene_symbol,
        "gene_disease_relations": gene_disease_relations,
        "gene_disease_summary": _get_gene_disease_summary(gene_disease_relations),
    }


def _get_gene_disease_summary(gene_disease_relations):
    """
        want it to look like 'Gene / Disease - 5 x autosomal dominant (supporting-definitive)'
    """

    moi_classification = defaultdict(Counter)
    for gdr in gene_disease_relations:
        for src in gdr.extra.get("sources", []):
            moi_classification[src["mode_of_inheritance"]][src["gencc_classification"]] += 1

    moi_summary = []
    for moi, classification_counts in moi_classification.items():
        classification_summary = ", ".join([f"{k}: {v}" for k, v in classification_counts.items()])
        moi_summary.append(f"{moi} - {classification_summary}")

    gene_disease_summary = "Gene / Disease: "
    if moi_summary:
        gene_disease_summary += " ".join(moi_summary)
    else:
        gene_disease_summary += "None"
    return gene_disease_summary
