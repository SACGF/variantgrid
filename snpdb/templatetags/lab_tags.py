from typing import Optional, Union

from django.conf import settings
from django.template.library import Library

from library.utils import group_data
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
def lab_picker(context, view_name: str, selected_lab: Optional[Union[Lab, int]] = None, all_option=False):
    labs_qs = Lab.valid_labs_qs(context.request.user, admin_check=True).select_related('organization')
    if isinstance(selected_lab, Lab):
        selected_lab = selected_lab.pk

    internal = labs_qs.filter(external=False)
    external = labs_qs.filter(external=True)

    return {
        "all_option": all_option,
        "internal": sorted(group_data(internal, lambda lab: (lab.organization, lab))),
        "external": sorted(group_data(external, lambda lab: (lab.organization, lab))),
        "lab_count": labs_qs.count(),
        "view_name": view_name,
        "selected_lab": selected_lab
    }
