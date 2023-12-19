from django.template import Library

from genes.models import GnomADGeneConstraint

register = Library()


@register.inclusion_tag("genes/tags/gnomad_gene_constraint_tag.html")
def gnomad_gene_constraint(gene_constraint: GnomADGeneConstraint):
    return {
        "gene_constraint": gene_constraint,
    }
