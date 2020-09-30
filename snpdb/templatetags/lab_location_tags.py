from django.template.library import Library

from snpdb.models import Lab

register = Library()


@register.inclusion_tag("snpdb/tags/lab_locations.html")
def lab_locations(labs=None, samples_only=False, involved_only=True,
                  center_lat=19.434403, center_long=37.238392, zoom_level=1):

    if labs is None:
        labs = Lab.objects.all()
    elif isinstance(labs, list):
        labs = Lab.objects.filter(pk__in=[lab.pk for lab in labs])

    labs = labs.filter(lat__isnull=False, long__isnull=False)
    if samples_only:
        labs = labs.filter(labproject__samplelabproject__sample__isnull=False)

    if involved_only:
        labs = labs.filter(labproject__involved=True)

    return {"labs": labs,
            "center_lat": center_lat,
            "center_long": center_long,
            "zoom_level": zoom_level}


@register.inclusion_tag("snpdb/tags/lab_families.html")
def lab_families(samples_only=False, involved_only=False, zoom_level=1):
    lab_info = Lab.objects.filter(labproject__families__gt=0)  # With families

    if samples_only:
        lab_info = lab_info.filter(labproject__samplelabproject__sample__isnull=False)

    if involved_only:
        lab_info = lab_info.filter(labproject__involved=True)

    return {"lab_info": lab_info,
            "zoom_level": zoom_level}
