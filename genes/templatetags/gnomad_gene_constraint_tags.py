from django.template import Library

from genes.models import GeneSymbol, GnomADGeneConstraint

register = Library()


@register.inclusion_tag("genes/tags/gnomad_gene_constraint_tag.html")
def gnomad_gene_constraint(gene_symbol: GeneSymbol):
    num_transcripts = GnomADGeneConstraint.objects.filter(gene_symbol=gene_symbol).count()
    gene_constraint = GnomADGeneConstraint.objects.filter(gene_symbol=gene_symbol).order_by("-mane_select").first()

    return {
        "gene_symbol": gene_symbol,
        "gene_constraint": gene_constraint,
        "num_transcripts": num_transcripts,
    }
