"""
The Effect node's rules over UTRAnnotator / SpliceVault text (annotation/effect_filters.py), run as the Postgres
regexes they are. Annotation strings are as UTRAnnotator writes them on GRCh38 (#2561).
"""
from django.db import connection
from django.test import TestCase

from annotation.effect_filters import (
    SPLICEVAULT_ALL_OUT_OF_FRAME_REGEX,
    SPLICEVAULT_ANY_OUT_OF_FRAME_REGEX,
    UTR_HIGH_IMPACT_RULES,
)

UAUG_GAIN = "5_prime_UTR_premature_start_codon_gain_variant"
UAUG_LOSS = "5_prime_UTR_premature_start_codon_loss_variant"
USTOP_LOSS = "5_prime_UTR_uORF_stop_codon_loss_variant"
UFRAMESHIFT = "5_prime_UTR_uORF_frameshift_variant"
USTOP_GAIN = "5_prime_UTR_uORF_stop_codon_gain_variant"


def _regex_matches(value: str, regex: str, case_insensitive: bool = False) -> bool:
    operator = "~*" if case_insensitive else "~"
    with connection.cursor() as cursor:
        cursor.execute(f"SELECT %s {operator} %s", [value, regex])
        return cursor.fetchone()[0]


class UTRHighImpactTest(TestCase):
    @staticmethod
    def _high_impact(consequence: str, annotation: str) -> bool:
        """ utr_high_impact_q, one rule per consequence """
        return any(consequence_term in consequence and _regex_matches(annotation, regex, case_insensitive=True)
                   for consequence_term, regex in UTR_HIGH_IMPACT_RULES)

    def test_uaug_gain_needs_oorf_and_strong_or_moderate_kozak(self):
        prefix = "CapDistanceToStart=111:DistanceToCDS=4:DistanceToStop=126:KozakContext=GTCATGG:"
        self.assertTrue(self._high_impact(UAUG_GAIN, prefix + "KozakStrength=Strong:type=OutOfFrame_oORF"))
        self.assertTrue(self._high_impact(UAUG_GAIN, prefix + "KozakStrength=Moderate:type=inFrame_oORF"))
        self.assertFalse(self._high_impact(UAUG_GAIN, prefix + "KozakStrength=Weak:type=OutOfFrame_oORF"))
        self.assertFalse(self._high_impact(UAUG_GAIN, prefix + "KozakStrength=Strong:type=uORF"))
        self.assertFalse(self._high_impact(UAUG_GAIN, "KozakContext=-:KozakStrength=:type=OutOfFrame_oORF"))

    def test_uaug_loss_needs_strong_kozak(self):
        prefix = "CapDistanceToStart=42:DistanceToCDS=89:DistanceToStop=21:Evidence=False:KozakContext=GCCATGG:"
        self.assertTrue(self._high_impact(UAUG_LOSS, prefix + "KozakStrength=Strong:type=uORF"))
        self.assertFalse(self._high_impact(UAUG_LOSS, prefix + "KozakStrength=Moderate:type=OutOfFrame_oORF"))

    def test_ustop_loss_needs_no_alternative_stop(self):
        suffix = ":Evidence=False:FrameWithCDS=outOfFrame:KozakContext=CGGATGG:KozakStrength=Moderate"
        self.assertTrue(self._high_impact(USTOP_LOSS, "AltStop=False:AltStopDistanceToCDS=NA" + suffix))
        self.assertFalse(self._high_impact(USTOP_LOSS, "AltStop=True:AltStopDistanceToCDS=113" + suffix))

    def test_uframeshift_needs_oorf(self):
        record = ("Evidence=True:KozakContext=GAGATGG:KozakStrength=Strong:alt_type={alt_type}:alt_type_length=270:"
                  "ref_StartDistanceToCDS=153:ref_type=OutOfFrame_oORF:ref_type_length=114")
        self.assertTrue(self._high_impact(UFRAMESHIFT, record.format(alt_type="OutOfFrame_oORF")))
        # ref_type is the uORF before the variant - only what the frameshift creates counts
        self.assertFalse(self._high_impact(UFRAMESHIFT, record.format(alt_type="uORF")))

    def test_ustop_gain_is_not_high_impact(self):
        annotation = ("Evidence=True:KozakContext=AGTATGA:KozakStrength=Strong:newSTOPDistanceToCDS=133:"
                      "ref_StartDistanceToCDS=175:ref_type=OutOfFrame_oORF")
        self.assertFalse(self._high_impact(USTOP_GAIN, annotation))

    def test_conditions_hold_within_one_effect(self):
        """ A weak-Kozak oORF beside a strong-Kozak uORF is two effects, neither high impact """
        weak_oorf = "CapDistanceToStart=68:DistanceToCDS=257:KozakStrength=Weak:type=OutOfFrame_oORF"
        strong_uorf = "CapDistanceToStart=68:DistanceToCDS=259:KozakStrength=Strong:type=uORF"
        self.assertFalse(self._high_impact(UAUG_GAIN, f"{weak_oorf}&{strong_uorf}"))
        strong_oorf = strong_uorf.replace("type=uORF", "type=OutOfFrame_oORF")
        self.assertTrue(self._high_impact(UAUG_GAIN, f"{weak_oorf}&{strong_oorf}"))

    def test_records_matched_to_their_consequence(self):
        """ A uAUG gain and loss on one variant: the loss record (it has Evidence) isn't read as a gain """
        consequence = f"{UAUG_GAIN}&{UAUG_LOSS}"
        weak_gain = "DistanceToCDS=257:KozakStrength=Weak:type=uORF"
        moderate_oorf_loss = "DistanceToCDS=259:Evidence=False:KozakStrength=Moderate:type=OutOfFrame_oORF"
        self.assertFalse(self._high_impact(consequence, f"{weak_gain}&{moderate_oorf_loss}"))


class SpliceVaultOutOfFrameTest(TestCase):
    def test_all_top_events_out_of_frame(self):
        self.assertTrue(_regex_matches("Frameshift:4/4", SPLICEVAULT_ALL_OUT_OF_FRAME_REGEX))
        self.assertTrue(_regex_matches("Frameshift:1/1", SPLICEVAULT_ALL_OUT_OF_FRAME_REGEX))
        self.assertFalse(_regex_matches("Frameshift:3/4", SPLICEVAULT_ALL_OUT_OF_FRAME_REGEX))
        self.assertFalse(_regex_matches("Frameshift:1/10", SPLICEVAULT_ALL_OUT_OF_FRAME_REGEX))
        self.assertFalse(_regex_matches("Frameshift:0/0", SPLICEVAULT_ALL_OUT_OF_FRAME_REGEX))

    def test_any_top_event_out_of_frame(self):
        self.assertTrue(_regex_matches("Frameshift:1/4", SPLICEVAULT_ANY_OUT_OF_FRAME_REGEX))
        self.assertFalse(_regex_matches("Frameshift:0/4", SPLICEVAULT_ANY_OUT_OF_FRAME_REGEX))
