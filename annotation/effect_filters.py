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

# 'Frameshift:N/M' - N of SpliceVault's M top events are out of frame. ALL needs N = M (and N > 0)
SPLICEVAULT_ALL_OUT_OF_FRAME_REGEX = r"^Frameshift:([1-9][0-9]*)/\1$"
SPLICEVAULT_ANY_OUT_OF_FRAME_REGEX = r"^Frameshift:[1-9]"


def splicevault_out_of_frame_q(all_events: bool, prefix: str = "variantannotation__") -> Q:
    regex = SPLICEVAULT_ALL_OUT_OF_FRAME_REGEX if all_events else SPLICEVAULT_ANY_OUT_OF_FRAME_REGEX
    return Q(**{f"{prefix}splicevault_out_of_frame_events__regex": regex})


# utr_annotator_annotation holds one record per 5' UTR effect, '&' separated, each 'key=value' pairs joined by ':'
# in no fixed order (annotation/vep_plugin_display.py:utr_annotator_effects). A record doesn't say which consequence
# it is for, but each consequence writes its own keys: a uAUG gain has 'type' and no 'Evidence', a uAUG loss 'type'
# and 'Evidence', a uSTOP loss 'AltStop', a uFrameshift 'alt_type'. Each lookahead stays inside one record ([^&]*),
# so every condition has to hold for the same effect. Matched case-insensitively: 'inFrame_oORF' and 'InFrame_oORF'
# both occur.
_RECORD_START = "(^|&)"
_OORF = "[a-z]*_oorf"  # an overlapping ORF in either frame - Martin-Geary 2025 doesn't restrict the frame
_STRONG_OR_MODERATE = "strong|moderate"


def _has_pair(key: str, value_regex: str) -> str:
    return f"(?=([^&]*:)?{key}=({value_regex})(:|&|$))"


def _lacks_key(key: str) -> str:
    return f"(?!([^&]*:)?{key}=)"


# (consequence, record regex) - Martin-Geary 2025's high-impact 5' UTR variants
UTR_HIGH_IMPACT_RULES: tuple[tuple[str, str], ...] = (
    # uAUG gain creating an oORF, strong or moderate Kozak
    ("5_prime_UTR_premature_start_codon_gain_variant",
     _RECORD_START + _has_pair("type", _OORF) + _has_pair("KozakStrength", _STRONG_OR_MODERATE)
     + _lacks_key("Evidence")),
    # uSTOP loss with no alternative stop before the CDS - so the uORF now overlaps it - strong or moderate Kozak
    ("5_prime_UTR_uORF_stop_codon_loss_variant",
     _RECORD_START + _has_pair("AltStop", "false") + _has_pair("KozakStrength", _STRONG_OR_MODERATE)),
    # uAUG loss, strong Kozak
    ("5_prime_UTR_premature_start_codon_loss_variant",
     _RECORD_START + _has_pair("type", "[^:&]*") + _has_pair("Evidence", "[^:&]*")
     + _has_pair("KozakStrength", "strong")),
    # uFrameshift creating an oORF, strong or moderate Kozak
    ("5_prime_UTR_uORF_frameshift_variant",
     _RECORD_START + _has_pair("alt_type", _OORF) + _has_pair("KozakStrength", _STRONG_OR_MODERATE)),
)


def utr_high_impact_q(prefix: str = "variantannotation__") -> Q:
    return reduce(operator.or_, (
        Q(**{f"{prefix}utr_annotator_consequence__contains": consequence,
             f"{prefix}utr_annotator_annotation__iregex": record_regex})
        for consequence, record_regex in UTR_HIGH_IMPACT_RULES
    ))
