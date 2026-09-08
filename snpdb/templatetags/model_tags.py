from typing import Optional

from django.template.library import Library
from django.utils.html import escape
from django.utils.safestring import mark_safe

from patients.models_enums import Sex
from snpdb.models import DETECTED_SEX_HELP, CohortGenotypeStats, Duo, Quad, Trio

register = Library()


@register.inclusion_tag("snpdb/tags/trio_table.html", takes_context=False)
def trio_table(trio: Trio):
    """ Trio details + a row per family member - the proband is affected by definition """
    members = [
        {"role": "Proband", "cohort_sample": trio.proband, "affected": True},
        {"role": "Mother", "cohort_sample": trio.mother, "affected": trio.mother_affected},
        {"role": "Father", "cohort_sample": trio.father, "affected": trio.father_affected},
    ]
    return {"trio": trio, "members": members}


@register.inclusion_tag("snpdb/tags/quad_table.html", takes_context=False)
def quad_table(quad: Quad):
    """ Quad details + a row per family member - the proband is affected by definition """
    members = [
        {"role": "Proband", "cohort_sample": quad.proband, "affected": True},
        {"role": "Mother", "cohort_sample": quad.mother, "affected": quad.mother_affected},
        {"role": "Father", "cohort_sample": quad.father, "affected": quad.father_affected},
        {"role": "Sibling", "cohort_sample": quad.sibling, "affected": quad.sibling_affected},
    ]
    return {"quad": quad, "members": members}


@register.inclusion_tag("snpdb/tags/duo_table.html", takes_context=False)
def duo_table(duo: Duo):
    """ Duo details + a row per family member - the proband is affected by definition """
    members = [
        {"role": "Proband", "cohort_sample": duo.proband, "affected": True},
        {"role": duo.relationship_label, "cohort_sample": duo.parent, "affected": duo.parent_affected},
    ]
    return {"duo": duo, "members": members}


@register.simple_tag
def trio_short_description(trio: Trio):
    params = (escape(trio.mother_details), escape(trio.father_details), escape(trio.proband))
    return mark_safe("<b>M:</b> %s/<b>F:</b> %s/<b>P:</b> %s" % params)


@register.simple_tag
def detected_sex_help(stats: Optional[CohortGenotypeStats] = None):
    """ Tooltip explaining how Sample.detected_sex is worked out - hang it off a label or column header.
        Pass the sample's stats to have the chrX counts it was called on spelled out too """
    if stats:
        return f"{DETECTED_SEX_HELP} {stats.chrx_sex_detail}"
    return DETECTED_SEX_HELP


@register.inclusion_tag("snpdb/tags/detected_sex.html", takes_context=False)
def detected_sex(stats: Optional[CohortGenotypeStats]):
    """ Detected sex as a table cell value, with this sample's chrX counts on hover """
    sex = stats.chrx_sex_guess if stats else Sex.UNKNOWN
    detail = stats.chrx_sex_detail if stats else "No chrX genotype counts for this sample"
    return {"sex": sex, "detail": detail}
