"""
Template filters that describe a Variant for a reader. Entry point: sv_size, the friendly size of a
symbolic structural variant ("11 Mb deletion") - in words, what the grids' badge shows
(VariantGridFormat.variantKind in variantgrid_formats.js).
"""
from django import template

from library.genomics import format_bp
from library.genomics.vcf_enums import VCFSymbolicAllele
from snpdb.models.models_variant import Variant

register = template.Library()

_SV_KIND_NAMES = {
    VCFSymbolicAllele.CNV: "CNV",
    VCFSymbolicAllele.DEL: "deletion",
    VCFSymbolicAllele.DUP: "duplication",
    VCFSymbolicAllele.INS: "insertion",
    VCFSymbolicAllele.INV: "inversion",
}


@register.filter
def sv_size(variant: Variant) -> str:
    """ '11 Mb deletion' for a <DEL> with svlen -11000000, '' for anything without a symbolic SV alt """
    kind = _SV_KIND_NAMES.get(variant.alt.seq)
    if not (kind and variant.svlen):
        return ""
    return f"{format_bp(abs(variant.svlen))} {kind}"
