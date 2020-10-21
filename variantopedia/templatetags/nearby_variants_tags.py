from django.template import Library

from annotation.models import AnnotationVersion
from snpdb.models import Variant
from variantopedia.interesting_nearby import get_nearby_summaries

register = Library()


@register.inclusion_tag("variantopedia/tags/nearby_variants.html", takes_context=True)
def nearby_variants(context, variant: Variant, annotation_version: AnnotationVersion,
                    clinical_significance: bool = True, distance: int = 50):
    user = context["user"]
    context.update({
        'variant': variant,
        "distance": distance,
    })
    context.update(get_nearby_summaries(user, variant, annotation_version,
                                        distance=distance, clinical_significance=clinical_significance))
    return context
