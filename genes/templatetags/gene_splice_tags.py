from typing import Optional

from django.template import Library

from genes.gene_splice import SpliceEventVariant
from snpdb.models import GenomeBuild

register = Library()


@register.inclusion_tag("genes/tags/splice_event.html")
def splice_event(splice_event_variant: SpliceEventVariant, genome_build: Optional[GenomeBuild] = None):
    """ The gene and the junction's name, in place of the storage coordinate a splice Variant
        formats as - @see snpdb.gene_level_variants. Given a build, the junction's breakpoints in it
        also get an IGV link (@see renderIgvLocusLinks in igv.js), and a junction CIViC records links
        to it (#1909) """
    igv_locus = splice_event_variant.junction_locus(genome_build) if genome_build else None
    splice_event = splice_event_variant.splice_event
    return {
        "splice_event_variant": splice_event_variant,
        "igv_locus": igv_locus,
        "civic_url": splice_event.civic_url if splice_event else None,
    }
