"""
Chips for the kinds of call a sample or VCF holds (Small, CNV, SV, Fusion, Gene CN, Splice - snpdb/variant_contents.py).

`variant_contents` draws one sample's or VCF's contents, as the grids' TableFormat.variantContents does.
`sample_variant_contents` / `vcf_variant_contents` look the contents up for a template's own list in one query, to
be read per row with `|get_item`. `variant_contents_help` is the legend, for a label's or column header's help.
"""
from collections.abc import Iterable
from typing import Optional

from django.template.library import Library
from django.utils.html import format_html_join
from django.utils.safestring import SafeString, mark_safe

from snpdb.models.models_vcf import VCF, Sample
from snpdb.variant_contents import (
    VARIANT_CONTENTS_HELP_LINES,
    VariantContents,
    get_sample_variant_contents,
    get_vcf_variant_contents,
    variant_contents_chips,
)

register = Library()


@register.inclusion_tag("snpdb/tags/variant_contents.html")
def variant_contents(contents: Optional[VariantContents], counts: bool = True, show_counts: bool = False,
                     placeholder: bool = False):
    """ counts=False for a VCF, whose numbers are sums over its samples. show_counts writes the count on the
        chip rather than only its tooltip. placeholder says why there are no chips (not counted yet, or
        counted and empty), for a page with a row to fill - a table cell is left empty """
    return {
        "chips": variant_contents_chips(contents, counts=counts),
        "show_counts": counts and show_counts,
        "placeholder": placeholder,
        "counted": contents is not None,
    }


@register.simple_tag
def sample_variant_contents(samples: Iterable[Sample]) -> dict[int, VariantContents]:
    return get_sample_variant_contents(samples)


@register.simple_tag
def vcf_variant_contents(vcfs: Iterable[VCF]) -> dict[int, VariantContents]:
    return get_vcf_variant_contents(vcfs)


@register.simple_tag
def variant_contents_help() -> SafeString:
    return format_html_join(mark_safe("<br>"), "{}", ((line,) for line in VARIANT_CONTENTS_HELP_LINES))
