from django.db.models import Count

from annotation.models import AnnotationRangeLock


def check_annotation_status() -> dict:
    # I am going to make this a warning for a while - before making it an error
    ARL_DUPE_ERROR = False
    # see https://github.com/SACGF/variantgrid_shariant/issues/177
    arl_qs = AnnotationRangeLock.objects.values('version_id', 'min_variant_id', 'max_variant_id')
    arl_qs = arl_qs.annotate(id_count=Count('id')).filter(id_count__gt=1)
    arl_message = "Remove duplicate AnnotationRangeLock objects"
    has_dupes = arl_qs.exists()
    if ARL_DUPE_ERROR:
        arl_data = {
            "valid": has_dupes,
            "fix": arl_message,
        }
    else:
        arl_data = {
            "valid": True,  # Just a warning
            "warning": arl_message + " (this is a just a warning but will become an error in the future)"
        }
    annotation_status = {
        "AnnotationRangeLock duplicates": arl_data,
    }
    return annotation_status
