from django.template import Library

from genes.gene_splice import SpliceEventVariant

register = Library()


@register.inclusion_tag("genes/tags/splice_event.html")
def splice_event(splice_event_variant: SpliceEventVariant):
    """ The gene and the junction's name, in place of the storage coordinate a splice Variant
        formats as - @see snpdb.gene_level_variants """
    return {"splice_event_variant": splice_event_variant}
