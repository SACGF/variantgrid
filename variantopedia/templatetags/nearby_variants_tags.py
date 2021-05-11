from django.conf import settings
from django.template import Library

from annotation.models import AnnotationVersion, GenomeBuild
from snpdb.models import Variant
from variantopedia.interesting_nearby import get_nearby_summaries, variant_interesting_summary

register = Library()


@register.inclusion_tag("variantopedia/tags/nearby_variants_tag.html", takes_context=True)
def nearby_variants(context, variant: Variant, annotation_version: AnnotationVersion,
                    clinical_significance: bool = True):
    distance: int = settings.VARIANT_DETAILS_NEARBY_RANGE
    user = context["user"]
    context.update({
        'variant': variant,
        "annotation_version": annotation_version,
        "distance": distance,
        "distance_str": str(distance)
    })
    context.update(get_nearby_summaries(user, variant, annotation_version,
                                        distance=distance, clinical_significance=clinical_significance))
    return context


@register.simple_tag(takes_context=True)
def nearby_summary(context, obj, genome_build: GenomeBuild):
    """ Only does anything for Variants """
    user = context["user"]
    summary = ""
    if isinstance(obj, Variant):
        summary = variant_interesting_summary(user, obj, genome_build)
    return summary
