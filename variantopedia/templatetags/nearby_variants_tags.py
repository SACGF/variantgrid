from django.template import Library

from annotation.models import AnnotationVersion
from snpdb.models import Variant

register = Library()


@register.inclusion_tag("variantopedia/tags/nearby_variants_tag.html", takes_context=True)
def nearby_variants(context, variant: Variant, annotation_version: AnnotationVersion,
                    clinical_significance: bool = True):
    context.update({
        'variant': variant,
        "annotation_version": annotation_version,
    })
    return context
