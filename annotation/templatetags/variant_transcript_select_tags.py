from django.template import Library

from genes.hgvs.hgvs import abbreviate_hgvs_sequences
from genes.hgvs.phgvs_canonical import abbreviate_p_hgvs_sequences

register = Library()


@register.inclusion_tag("annotation/tags/variant_transcript_select.html")
def variant_transcript_select(transcript_select_jquery, vts,
                              transcript_func=None,
                              gene_func=None,
                              show_all_transcript_details=False,
                              empty_transcript_option=None):
    return {
        "transcript_select_jquery": transcript_select_jquery,
        "transcript_func": transcript_func,
        "gene_func": gene_func,
        "vts": vts,
        "variant_transcript_annotations": vts.variant_transcript_annotations_dict,
        "show_all_transcript_details": show_all_transcript_details,
        "empty_transcript_option": empty_transcript_option
    }


@register.filter
def abbreviate_hgvs(hgvs):
    return abbreviate_hgvs_sequences(hgvs)


@register.filter
def abbreviate_p_hgvs(p_hgvs):
    return abbreviate_p_hgvs_sequences(p_hgvs)
