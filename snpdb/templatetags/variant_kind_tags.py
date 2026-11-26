"""
The variant-kind badge (DEL 20 kb, DUP, LOSS, fusion...) beside a symbolic variant a template names.

`variant_kind_badge` emits the placeholder that `VariantGridFormat.variantKindBadge` (variantgrid_formats.js)
replaces, so templates show the same badge as the grids. Empty for a small variant. `whole_chromosome_iscn` is
the label a grid puts beside a whole chromosome's badge ("+21").
"""
from typing import Any

from django.template.library import Library
from django.utils.html import format_html, format_html_join

from library.genomics.vcf_enums import VCFSymbolicAllele
from snpdb.models.models_variant import Allele, Variant

register = Library()


@register.simple_tag
def variant_kind_badge(obj: Any) -> str:
    """ A Variant, or an Allele - whose kind and length are the same in every build, so its badge's
        tooltip leaves out the coordinate. Anything else gets no badge """
    if isinstance(obj, Allele):
        variant = obj.variants.select_related("alt").first()
        coordinate = {}
    elif isinstance(obj, Variant):
        variant = obj
        coordinate = {"chrom": variant.locus.contig.name, "position": variant.locus.position}
    else:
        return ""

    if not (variant and variant.alt.is_symbolic):
        return ""
    attributes = {"alt": variant.alt.seq, "svlen": variant.svlen or 0, **coordinate}
    return format_html("<span{}></span>",
                       format_html_join("", ' data-variant-kind-{}="{}"', attributes.items()))


@register.simple_tag
def whole_chromosome_iscn(variant: Variant) -> str:
    """ A whole chromosome's ISCN ("+21") in its badge's colour, as the grid labels it
        (VariantGridFormat.representativeVariantLabel). Empty for any other variant """
    if not (iscn := variant.whole_chromosome_iscn):
        return ""
    css = "rv-chrom-gain" if variant.alt.seq == VCFSymbolicAllele.DUP_CHROM else "rv-chrom-loss"
    return format_html("<span class='rv-gene {}'>{}</span>", css, iscn)
