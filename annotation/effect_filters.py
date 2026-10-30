"""
Q filters for the Effect (Damage) node over the categorical plugin annotations, where "damaging" is a rule over
the text VEP wrote rather than a score threshold (#2561): SpliceVault's out-of-frame event count and the
UTRAnnotator high-impact 5' UTR rule. The rules and their sources are listed on the Pathogenicity Thresholds page
(`annotation/pathogenicity_predictions.py:EFFECT_NODE_FILTERS`).

Entry points: splicevault_out_of_frame_q, utr_high_impact_q - each takes the path prefix to VariantAnnotation.
"""
import operator
from functools import reduce

from django.db.models import Q

# splicevault_out_of_frame_events is 'Frameshift:N/M' - N of SpliceVault's M top events are out of frame. All of
# them needs N = M (and N > 0), Dawes 2023's "consistent with null outcomes" case
SPLICEVAULT_ALL_OUT_OF_FRAME_REGEX = r"^Frameshift:([1-9][0-9]*)/\1$"
# splicevault_top_events is '&'-joined 'TopN:<event>:<position>:<% of samples>:<frame>', most frequent first. The
# Top1 event being out of frame weights by frequency, which the N/M count can't: a 12% frameshift ahead of three
# rare in-frame events is only 1/4, but the variant's usual outcome is still a frameshift. "Any event out of frame"
# isn't offered - with 4 top events at a site SpliceAI already predicts lost, nearly every variant has one
SPLICEVAULT_TOP_EVENT_OUT_OF_FRAME_REGEX = r"^Top1:[^&]*:Frameshift(&|$)"


def splicevault_out_of_frame_q(top_event_only: bool, prefix: str = "variantannotation__") -> Q:
    """ top_event_only: SpliceVault's most frequent event is out of frame, else all of its top events are """
    if top_event_only:
        return Q(**{f"{prefix}splicevault_top_events__regex": SPLICEVAULT_TOP_EVENT_OUT_OF_FRAME_REGEX})
    return Q(**{f"{prefix}splicevault_out_of_frame_events__regex": SPLICEVAULT_ALL_OUT_OF_FRAME_REGEX})


# utr_annotator_annotation holds one record per uORF effect, each 'key=value' pairs joined by ':'
# (annotation/vep_plugin_display.py:utr_annotator_effects). The plugin writes the effects in a fixed order (uAUG gain,
# uSTOP loss, uAUG loss, uSTOP gain, uFrameshift - the order of utr_annotator_consequence), records of one effect
# joined by ',' and effects by '&', but VEP's VCF output turns ',' into '&' as well, so the stored string no longer
# says where one effect's records end. Each effect writes its own set of keys, so a record says which effect it is:
# a uAUG gain has 'type' and never 'Evidence', a uAUG loss 'type' and always 'Evidence' ('NA' without an evidence
# file), a uSTOP loss 'AltStop', a uFrameshift 'alt_type', a uSTOP gain 'newSTOPDistanceToCDS'. Each lookahead stays
# inside one record ([^&]*), so every condition has to hold for the same effect. Matched case-insensitively:
# 'inFrame_oORF' and 'InFrame_oORF' both occur.
_RECORD_START = "(^|&)"
_STRONG_OR_MODERATE = "strong|moderate"


def _has_pair(key: str, value_regex: str) -> str:
    return f"(?=([^&]*:)?{key}=({value_regex})(:|&|$))"


def _lacks_key(key: str) -> str:
    return f"(?!([^&]*:)?{key}=)"


def utr_high_impact_rules(out_of_frame_only: bool) -> tuple[tuple[str, str], ...]:
    """ (consequence, record regex) for Martin-Geary 2025's high-impact 5' UTR variants. Martin-Geary doesn't
        restrict the oORF's frame; out_of_frame_only drops an oORF in frame with the CDS (an N-terminal extension,
        not necessarily loss of function - Ellingford 2022) from the three effects that create one """
    oorf = "outofframe_oorf" if out_of_frame_only else "[a-z]*_oorf"
    ustop_lost_frame = _has_pair("FrameWithCDS", "outofframe") if out_of_frame_only else ""
    return (
        # uAUG gain creating an oORF, strong or moderate Kozak
        ("5_prime_UTR_premature_start_codon_gain_variant",
         _RECORD_START + _has_pair("type", oorf) + _has_pair("KozakStrength", _STRONG_OR_MODERATE)
         + _lacks_key("Evidence")),
        # uSTOP loss with no alternative stop before the CDS - so the uORF now overlaps it - strong or moderate Kozak
        ("5_prime_UTR_uORF_stop_codon_loss_variant",
         _RECORD_START + _has_pair("AltStop", "false") + ustop_lost_frame
         + _has_pair("KozakStrength", _STRONG_OR_MODERATE)),
        # uAUG loss, strong Kozak - whatever the frame
        ("5_prime_UTR_premature_start_codon_loss_variant",
         _RECORD_START + _has_pair("type", "[^:&]*") + _has_pair("Evidence", "[^:&]*")
         + _has_pair("KozakStrength", "strong")),
        # uFrameshift creating an oORF, strong or moderate Kozak
        ("5_prime_UTR_uORF_frameshift_variant",
         _RECORD_START + _has_pair("alt_type", oorf) + _has_pair("KozakStrength", _STRONG_OR_MODERATE)),
    )


def utr_high_impact_q(out_of_frame_only: bool, prefix: str = "variantannotation__") -> Q:
    return reduce(operator.or_, (
        Q(**{f"{prefix}utr_annotator_consequence__contains": consequence,
             f"{prefix}utr_annotator_annotation__iregex": record_regex})
        for consequence, record_regex in utr_high_impact_rules(out_of_frame_only)
    ))
