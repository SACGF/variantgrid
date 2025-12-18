from django.template.library import Library

from snpdb.models import Lab
from variantgrid import settings

register = Library()


from django import template
from django.template.loader import render_to_string
from django.conf import settings

register = template.Library()


@register.simple_tag
def lab_locations(labs=None, samples_only=False, involved_only=True,
                  center_lat=19.434403, center_long=37.238392, zoom_level=1):
    # Made this a simple tag, so that if Leaflet is not loaded, we don't load any Leaflet tags

    if "leaflet" not in settings.INSTALLED_APPS:
        return ""

    if labs is None:
        labs = Lab.objects.all()
    elif isinstance(labs, list):
        labs = Lab.objects.filter(pk__in=[lab.pk for lab in labs])

    labs = labs.filter(lat__isnull=False, long__isnull=False)

    if samples_only:
        labs = labs.filter(labproject__samplelabproject__sample__isnull=False)
    if involved_only:
        labs = labs.filter(labproject__involved=True)

    if not labs.exists():
        return ""

    return render_to_string(
        "snpdb/tags/lab_locations.html",
        {
            "labs": labs,
            "center_lat": center_lat,
            "center_long": center_long,
            "zoom_level": zoom_level,
        }
    )


@register.inclusion_tag("snpdb/tags/lab_families.html")
def lab_families(samples_only=False, involved_only=False, zoom_level=1):
    lab_info = Lab.objects.filter(labproject__families__gt=0)  # With families

    if samples_only:
        lab_info = lab_info.filter(labproject__samplelabproject__sample__isnull=False)

    if involved_only:
        lab_info = lab_info.filter(labproject__involved=True)

    return {"lab_info": lab_info,
            "zoom_level": zoom_level}
