from django.template import Library

from genes.models import GeneFusion

register = Library()


@register.inclusion_tag("genes/tags/gene_fusion.html")
def gene_fusion(fusion: GeneFusion):
    """ The partners and their direction, in place of the storage coordinate a fusion Variant
        formats as - @see snpdb.gene_level_variants """
    return {"gene_fusion": fusion}
