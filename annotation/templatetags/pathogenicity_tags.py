import json

from django import template
from django.utils.safestring import mark_safe

from annotation.pathogenicity_predictions import TOOLS

register = template.Library()


@register.simple_tag
def pathogenicity_thresholds() -> str:
    """JSON object {field: [low, high]} or {field: [threshold]} keyed by raw_field,
    consumed by variant_details.html CODE_THRESHOLDS for colour banding."""
    bands: dict[str, list[float]] = {}
    for t in TOOLS:
        if not t.raw_field or t.raw_pathogenic_threshold is None:
            continue
        if t.raw_max_benign_threshold is None:
            bands[t.raw_field] = [t.raw_pathogenic_threshold]
        else:
            bands[t.raw_field] = [t.raw_max_benign_threshold, t.raw_pathogenic_threshold]
    # Embedded inside <script> in variant_details.html — JSON is safe to inject.
    return mark_safe(json.dumps(bands))
