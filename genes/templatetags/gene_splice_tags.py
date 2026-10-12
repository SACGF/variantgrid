from typing import Optional

from django.template import Library

from genes.gene_splice import SpliceEventVariant
from snpdb.models import GenomeBuild

register = Library()


@register.inclusion_tag("genes/tags/splice_event.html")
def splice_event(splice_event_variant: SpliceEventVariant, genome_build: Optional[GenomeBuild] = None):
    """ The gene and the junction's name, in place of the storage coordinate a splice Variant
        formats as - @see snpdb.gene_level_variants. Given a build, it also carries the junction's
        links, for a page with no Quick Links to put them in: IGV at its breakpoints in that build
        (@see renderIgvLocusLinks in igv.js), and CIViC for a junction CIViC records (#1909) """
    context = {"splice_event_variant": splice_event_variant}
    if genome_build:
        splice_event = splice_event_variant.splice_event
        context["igv_locus"] = splice_event_variant.junction_locus(genome_build)
        context["civic_url"] = splice_event.civic_url if splice_event else None
    return context
