from django.template import Library

register = Library()

@register.inclusion_tag("genes/tags/gnomad_gene_constraint.html")
def gnomad_gene_constraint(gene_symbol):
    return {"gene_symbol": gene_symbol}
