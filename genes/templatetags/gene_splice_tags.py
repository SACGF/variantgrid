from typing import Optional

from django.template import Library

from genes.gene_splice import SpliceEventVariant
from snpdb.models import GenomeBuild

register = Library()


@register.inclusion_tag("genes/tags/splice_event.html")
def splice_event(splice_event_variant: SpliceEventVariant, genome_build: Optional[GenomeBuild] = None):
    """ The gene and the junction's name, in place of the storage coordinate a splice Variant
        formats as - @see snpdb.gene_level_variants. Given a build, the junction's breakpoints in it
        also get an IGV link (@see renderIgvLocusLinks in grid.js) """
    igv_locus = splice_event_variant.junction_locus(genome_build) if genome_build else None
    return {
        "splice_event_variant": splice_event_variant,
        "igv_locus": igv_locus,
    }
