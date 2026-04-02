from django.template.library import Library
from django.utils.html import escape
from django.utils.safestring import mark_safe

from snpdb.models import Trio, Quad

register = Library()


@register.inclusion_tag("snpdb/tags/trio_table.html", takes_context=False)
def trio_table(trio: Trio):
    return {"trio": trio}


@register.inclusion_tag("snpdb/tags/quad_table.html", takes_context=False)
def quad_table(quad: Quad):
    return {"quad": quad}


@register.simple_tag
def trio_short_description(trio: Trio):
    params = (escape(trio.mother_details), escape(trio.father_details), escape(trio.proband))
    return mark_safe("<b>M:</b> %s/<b>F:</b> %s/<b>P:</b> %s" % params)
