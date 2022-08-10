from django.conf import settings
from django.template.library import Library

from snpdb.lab_picker import LabPickerData
from snpdb.models import Lab

register = Library()


@register.inclusion_tag("snpdb/tags/lab_card.html", takes_context=True)
def lab_card(context, lab: Lab, lab_link=True, org_link=True):
    user = context.request.user
    return {
        "user": user,
        "lab": lab,
        "org": lab.organization,
        "lab_link": lab_link,
        "org_link": org_link,
        "is_member": lab.is_member(user) or user.is_superuser,
        "shared_classifications": settings.VARIANT_CLASSIFICATION_STATS_USE_SHARED,
        "url_name_visible": context["url_name_visible"]
    }


@register.inclusion_tag("snpdb/tags/lab_picker.html", takes_context=True)
def lab_picker(context, data: LabPickerData):
    return {"data": data}
