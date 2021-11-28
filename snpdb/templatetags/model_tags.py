from django.template.library import Library
from django.utils.safestring import mark_safe

from snpdb.models import Trio

register = Library()


@register.inclusion_tag("snpdb/tags/trio_table.html", takes_context=False)
def trio_table(trio: Trio):
    return {"trio": trio}


@register.simple_tag
def trio_short_description(trio: Trio):
    params = (trio.mother_details, trio.father_details, trio.proband)
    return mark_safe("<b>M:</b> %s/<b>F:</b> %s/<b>P:</b> %s" % params)
