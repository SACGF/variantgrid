from typing import Optional

from django.template.library import Library
from django.utils.html import escape
from django.utils.safestring import mark_safe

from snpdb.models import DETECTED_SEX_HELP, CohortGenotypeStats, Duo, Quad, Sample, Trio
from snpdb.variant_contents import get_sample_variant_contents

register = Library()


def _with_variant_contents(members: list[dict]) -> list[dict]:
    """ Each family member's row gets the kinds of call its sample holds, looked up for the family at once """
    contents = get_sample_variant_contents([member["cohort_sample"].sample_id for member in members])
    for member in members:
        member["variant_contents"] = contents.get(member["cohort_sample"].sample_id)
    return members


@register.inclusion_tag("snpdb/tags/trio_table.html", takes_context=False)
def trio_table(trio: Trio):
    """ Trio details + a row per family member - the proband is affected by definition """
    members = [
        {"role": "Father", "cohort_sample": trio.father, "affected": trio.father_affected},
        {"role": "Mother", "cohort_sample": trio.mother, "affected": trio.mother_affected},
        {"role": "Proband", "cohort_sample": trio.proband, "affected": True},
    ]
    return {"trio": trio, "members": _with_variant_contents(members)}


@register.inclusion_tag("snpdb/tags/quad_table.html", takes_context=False)
def quad_table(quad: Quad):
    """ Quad details + a row per family member - the proband is affected by definition """
    members = [
        {"role": "Father", "cohort_sample": quad.father, "affected": quad.father_affected},
        {"role": "Mother", "cohort_sample": quad.mother, "affected": quad.mother_affected},
        {"role": "Proband", "cohort_sample": quad.proband, "affected": True},
        {"role": "Sibling", "cohort_sample": quad.sibling, "affected": quad.sibling_affected},
    ]
    return {"quad": quad, "members": _with_variant_contents(members)}


@register.inclusion_tag("snpdb/tags/duo_table.html", takes_context=False)
def duo_table(duo: Duo):
    """ Duo details + a row per family member - the proband is affected by definition """
    members = [
        {"role": "Proband", "cohort_sample": duo.proband, "affected": True},
        {"role": duo.relationship_label, "cohort_sample": duo.relative, "affected": duo.relative_affected},
    ]
    return {"duo": duo, "members": _with_variant_contents(members)}


@register.simple_tag
def trio_short_description(trio: Trio):
    params = (escape(trio.mother_details), escape(trio.father_details), escape(trio.proband))
    return mark_safe("<b>M:</b> %s/<b>F:</b> %s/<b>P:</b> %s" % params)


@register.simple_tag
def detected_sex_help(sample: Optional[Sample] = None, stats: Optional[CohortGenotypeStats] = None):
    """ Tooltip explaining how Sample.detected_sex is worked out - hang it off a label or column header.
        Pass the sample (and its stats, if already fetched) to have what it was detected from spelled out too """
    if sample:
        _, detail = sample.get_detected_sex(stats)
        return f"{DETECTED_SEX_HELP} {detail}"
    return DETECTED_SEX_HELP


@register.inclusion_tag("snpdb/tags/detected_sex.html", takes_context=False)
def detected_sex(sample: Sample, stats: Optional[CohortGenotypeStats] = None):
    """ Detected sex as a table cell value, with what it was detected from on hover. Pass stats already
        fetched for a table of samples to save the lookup """
    sex, detail = sample.get_detected_sex(stats)
    return {"sex": sex, "detail": detail}
