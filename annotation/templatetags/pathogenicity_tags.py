import json

from django import template
from django.conf import settings
from django.utils.safestring import mark_safe

from annotation.pathogenicity_predictions import COLOUR_BANDS, TOOLS

register = template.Library()


def _band(max_benign: float, pathogenic: float, magnitude: bool = False) -> dict:
    band = {"pathogenic": pathogenic}
    if max_benign is not None:
        band["benign"] = max_benign
    if magnitude:
        band["magnitude"] = True
    return band


@register.simple_tag
def pathogenicity_thresholds() -> str:
    """JSON object {field: {benign, pathogenic, magnitude}} keyed by annotation column, consumed by
    variant_details.html CODE_THRESHOLDS for colour banding. A band with no "benign" colours only the
    damaging end. Rankscores band off the site-wide settings rather than a per-tool cutoff."""
    bands: dict[str, dict] = {}
    for t in TOOLS:
        if t.rankscore_field:
            bands[t.rankscore_field] = _band(settings.ANNOTATION_MAX_BENIGN_RANKSCORE,
                                             settings.ANNOTATION_MIN_PATHOGENIC_RANKSCORE)
        if t.raw_field and t.raw_pathogenic_threshold is not None:
            bands[t.raw_field] = _band(t.raw_max_benign_threshold, t.raw_pathogenic_threshold)
    for cb in COLOUR_BANDS:
        if cb.pathogenic_threshold is None:
            continue
        for field in cb.fields:
            bands[field] = _band(cb.max_benign_threshold, cb.pathogenic_threshold, cb.magnitude)
    # Embedded inside <script> in variant_details.html — JSON is safe to inject.
    return mark_safe(json.dumps(bands))
