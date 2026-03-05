from django.template import Library

register = Library()


@register.inclusion_tag("genes/tags/hgnc_tag.html")
def hgnc_tag(hgnc, annotation_description):
    return {
        "hgnc": hgnc,
        "annotation_description": annotation_description,
    }
