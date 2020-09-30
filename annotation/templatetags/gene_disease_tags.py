from django.template import Library

register = Library()

@register.inclusion_tag("annotation/tags/gene_disease.html")
def gene_disease(gene_symbol):
    return {"gene_symbol": gene_symbol}
