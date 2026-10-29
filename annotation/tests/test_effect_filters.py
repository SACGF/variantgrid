"""
The Effect node's rules over UTRAnnotator / SpliceVault text (annotation/effect_filters.py), run as the Postgres
regexes they are. Annotation strings are as UTRAnnotator writes them on GRCh38 (#2561).
"""
from django.db import connection
from django.test import TestCase

from annotation.effect_filters import (
    SPLICEVAULT_ALL_OUT_OF_FRAME_REGEX,
    SPLICEVAULT_TOP_EVENT_OUT_OF_FRAME_REGEX,
    utr_high_impact_rules,
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


def _high_impact(consequence: str, annotation: str, out_of_frame_only: bool = False) -> bool:
    """ utr_high_impact_q, one rule per consequence """
    return any(consequence_term in consequence and _regex_matches(annotation, regex, case_insensitive=True)
               for consequence_term, regex in utr_high_impact_rules(out_of_frame_only))


class UTRHighImpactTest(TestCase):
    def test_uaug_gain_needs_oorf_and_strong_or_moderate_kozak(self):
        prefix = "CapDistanceToStart=111:DistanceToCDS=4:DistanceToStop=126:KozakContext=GTCATGG:"
        self.assertTrue(_high_impact(UAUG_GAIN, prefix + "KozakStrength=Strong:type=OutOfFrame_oORF"))
        self.assertTrue(_high_impact(UAUG_GAIN, prefix + "KozakStrength=Moderate:type=inFrame_oORF"))
        self.assertFalse(_high_impact(UAUG_GAIN, prefix + "KozakStrength=Weak:type=OutOfFrame_oORF"))
        self.assertFalse(_high_impact(UAUG_GAIN, prefix + "KozakStrength=Strong:type=uORF"))
        self.assertFalse(_high_impact(UAUG_GAIN, "KozakContext=-:KozakStrength=:type=OutOfFrame_oORF"))

    def test_uaug_loss_needs_strong_kozak(self):
        prefix = "CapDistanceToStart=42:DistanceToCDS=89:DistanceToStop=21:Evidence=False:KozakContext=GCCATGG:"
        self.assertTrue(_high_impact(UAUG_LOSS, prefix + "KozakStrength=Strong:type=uORF"))
        self.assertTrue(_high_impact(UAUG_LOSS, prefix + "KozakStrength=Strong:type=InFrame_oORF",
                                     out_of_frame_only=True))
        self.assertFalse(_high_impact(UAUG_LOSS, prefix + "KozakStrength=Moderate:type=OutOfFrame_oORF"))

    def test_ustop_loss_needs_no_alternative_stop(self):
        suffix = ":Evidence=False:FrameWithCDS=outOfFrame:KozakContext=CGGATGG:KozakStrength=Moderate"
        self.assertTrue(_high_impact(USTOP_LOSS, "AltStop=False:AltStopDistanceToCDS=NA" + suffix))
        self.assertFalse(_high_impact(USTOP_LOSS, "AltStop=True:AltStopDistanceToCDS=113" + suffix))

    def test_uframeshift_needs_oorf(self):
        record = ("Evidence=True:KozakContext=GAGATGG:KozakStrength=Strong:alt_type={alt_type}:alt_type_length=270:"
                  "ref_StartDistanceToCDS=153:ref_type=OutOfFrame_oORF:ref_type_length=114")
        self.assertTrue(_high_impact(UFRAMESHIFT, record.format(alt_type="OutOfFrame_oORF")))
        # ref_type is the uORF before the variant - only what the frameshift creates counts
        self.assertFalse(_high_impact(UFRAMESHIFT, record.format(alt_type="uORF")))

    def test_out_of_frame_only_drops_in_frame_oorfs(self):
        """ An oORF in frame with the CDS is an N-terminal extension - kept by "any frame", dropped otherwise """
        in_frame_gain = "CapDistanceToStart=141:DistanceToCDS=15:KozakStrength=Moderate:type=inFrame_oORF"
        in_frame_shift = "Evidence=False:KozakStrength=Moderate:alt_type=InFrame_oORF:ref_type=uORF"
        in_frame_stop_loss = "AltStop=False:AltStopDistanceToCDS=NA:Evidence=False:FrameWithCDS=inFrame:KozakStrength=Strong"
        for consequence, annotation in [(UAUG_GAIN, in_frame_gain), (UFRAMESHIFT, in_frame_shift),
                                        (USTOP_LOSS, in_frame_stop_loss)]:
            self.assertTrue(_high_impact(consequence, annotation))
            self.assertFalse(_high_impact(consequence, annotation, out_of_frame_only=True))

        out_of_frame_stop_loss = in_frame_stop_loss.replace("FrameWithCDS=inFrame", "FrameWithCDS=outOfFrame")
        self.assertTrue(_high_impact(USTOP_LOSS, out_of_frame_stop_loss, out_of_frame_only=True))
        out_of_frame_gain = in_frame_gain.replace("inFrame_oORF", "OutOfFrame_oORF")
        self.assertTrue(_high_impact(UAUG_GAIN, out_of_frame_gain, out_of_frame_only=True))

    def test_ustop_gain_is_not_high_impact(self):
        annotation = ("Evidence=True:KozakContext=AGTATGA:KozakStrength=Strong:newSTOPDistanceToCDS=133:"
                      "ref_StartDistanceToCDS=175:ref_type=OutOfFrame_oORF")
        self.assertFalse(_high_impact(USTOP_GAIN, annotation))

    def test_conditions_hold_within_one_effect(self):
        """ A weak-Kozak oORF beside a strong-Kozak uORF is two effects, neither high impact """
        weak_oorf = "CapDistanceToStart=68:DistanceToCDS=257:KozakStrength=Weak:type=OutOfFrame_oORF"
        strong_uorf = "CapDistanceToStart=68:DistanceToCDS=259:KozakStrength=Strong:type=uORF"
        self.assertFalse(_high_impact(UAUG_GAIN, f"{weak_oorf}&{strong_uorf}"))
        strong_oorf = strong_uorf.replace("type=uORF", "type=OutOfFrame_oORF")
        self.assertTrue(_high_impact(UAUG_GAIN, f"{weak_oorf}&{strong_oorf}"))

    def test_gain_and_loss_records_attributed_to_their_effect(self):
        """ A variant both gaining and losing a uAUG: the plugin writes the gain's records (never 'Evidence') then the
            loss's (always 'Evidence'), so a moderate-Kozak oORF the variant destroys isn't read as one it creates """
        consequence = f"{UAUG_GAIN}&{UAUG_LOSS}"
        gain = "CapDistanceToStart=54:DistanceToCDS=56:DistanceToStop=57:KozakContext=TCTATGG:KozakStrength={kozak}:type={type}"
        loss = ("CapDistanceToStart=54:DistanceToCDS=59:DistanceToStop=54:Evidence=False:KozakContext=TCTATGG:"
                "KozakStrength=Moderate:type=OutOfFrame_oORF")
        self.assertFalse(_high_impact(consequence, gain.format(kozak="Weak", type="OutOfFrame_oORF") + "&" + loss))
        self.assertFalse(_high_impact(consequence, gain.format(kozak="Moderate", type="uORF") + "&" + loss))
        self.assertTrue(_high_impact(consequence, gain.format(kozak="Moderate", type="OutOfFrame_oORF") + "&" + loss))


class SpliceVaultOutOfFrameTest(TestCase):
    def test_all_top_events_out_of_frame(self):
        self.assertTrue(_regex_matches("Frameshift:4/4", SPLICEVAULT_ALL_OUT_OF_FRAME_REGEX))
        self.assertTrue(_regex_matches("Frameshift:1/1", SPLICEVAULT_ALL_OUT_OF_FRAME_REGEX))
        self.assertFalse(_regex_matches("Frameshift:3/4", SPLICEVAULT_ALL_OUT_OF_FRAME_REGEX))
        self.assertFalse(_regex_matches("Frameshift:1/10", SPLICEVAULT_ALL_OUT_OF_FRAME_REGEX))
        self.assertFalse(_regex_matches("Frameshift:0/0", SPLICEVAULT_ALL_OUT_OF_FRAME_REGEX))

    def test_top_event_out_of_frame(self):
        """ The most frequent event decides, whatever the rarer ones are """
        top1_frameshift = ("Top1:CD:-32:12%:Frameshift&Top2:ES:2:7%:inFrame&Top3:CD:-96:0.4%:inFrame&"
                           "Top4:ES:2-3:0.3%:inFrame")
        self.assertTrue(_regex_matches(top1_frameshift, SPLICEVAULT_TOP_EVENT_OUT_OF_FRAME_REGEX))
        self.assertTrue(_regex_matches("Top1:ES:3-4:35%:Frameshift", SPLICEVAULT_TOP_EVENT_OUT_OF_FRAME_REGEX))
        top1_in_frame = ("Top1:ES:3:36%:inFrame&Top2:ES:3-4:34%:Frameshift&Top3:ES:2-3:0.3%:inFrame&"
                         "Top4:CA:+17:0.05%:Frameshift")
        self.assertFalse(_regex_matches(top1_in_frame, SPLICEVAULT_TOP_EVENT_OUT_OF_FRAME_REGEX))
