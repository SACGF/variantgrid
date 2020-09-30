from django.template import Library
import json

register = Library()


@register.inclusion_tag("snpdb/tags/model_fields_version_diff.html")
def model_fields_version_diff(model_dicts_by_version, versions):
    return {"model_dicts_by_version": model_dicts_by_version,
            "versions": versions}
