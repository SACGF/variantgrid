from django.template import Library

from genes.models import GeneCopyNumberEvent

register = Library()


@register.inclusion_tag("genes/tags/gene_copy_number.html")
def gene_copy_number(event: GeneCopyNumberEvent):
    """ The gene and which way it went, in place of the storage coordinate a copy number Variant
        formats as - @see snpdb.gene_level_variants """
    return {"gene_copy_number_event": event}
