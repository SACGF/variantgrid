from django.test import TestCase

from annotation.models.repeat_masker import (
    RepeatMaskerSummary, classify_repeat,
    SIMPLE_REPEAT, LOW_COMPLEXITY, SINE, LINE, LTR, DNA, RNA, OTHER,
)

# Real value from issue #1580 (NM_182961.4(SYNE1):c.5270_23300del), truncated
ISSUE_1580_VALUE = (
    "Tigger1&L1MA2&L1M1&AluSx3&L1M1&AluSz6&L1M1&AluY&L1M1&AluSx&L1M1&Tigger1&G-rich&"
    "Tigger1&MIRb&AluJr4&X7D_LINE&AluSx1&MIRb&L2a&(CTG)n&MIR&(T)n&A-rich&MER41C&"
    "MER34A&MER4D1&5S&AmnSINE1&(AC)n&MADE1&Charlie26a"
)


class RepeatMaskerTests(TestCase):
    def test_classify_repeat(self):
        cases = {
            "AluSx3": SINE, "AluY": SINE, "MIRb": SINE, "AmnSINE1": SINE,
            "L1M1": LINE, "L2a": LINE, "L3b": LINE, "HAL1b": LINE, "X7D_LINE": LINE,
            "(CTG)n": SIMPLE_REPEAT, "(T)n": SIMPLE_REPEAT,
            "A-rich": LOW_COMPLEXITY, "G-rich": LOW_COMPLEXITY, "GA-rich": LOW_COMPLEXITY,
            "LTR10B1": LTR, "THE1B": LTR, "MLT1A1": LTR,
            "MER41C": LTR, "MER34A": LTR,   # LTR-class MER
            "Tigger1": DNA, "MADE1": DNA, "Charlie26a": DNA,
            "MER30": DNA, "MER20": DNA,     # DNA-class MER
            "5S": RNA,
            "totally_unknown_thing": OTHER,
        }
        for name, expected in cases.items():
            self.assertEqual(classify_repeat(name), expected, name)

    def test_summary_empty(self):
        for value in (None, ""):
            summary = RepeatMaskerSummary.from_value(value)
            self.assertFalse(summary.is_multiple)
            self.assertEqual(summary.total, 0)
            self.assertEqual(summary.class_counts, [])

    def test_summary_single(self):
        summary = RepeatMaskerSummary.from_value("AluY")
        self.assertFalse(summary.is_multiple)
        self.assertEqual(summary.total, 1)
        self.assertEqual(summary.class_counts, [(SINE, 1)])

    def test_summary_issue_1580(self):
        summary = RepeatMaskerSummary.from_value(ISSUE_1580_VALUE)
        self.assertTrue(summary.is_multiple)
        self.assertEqual(summary.total, len(ISSUE_1580_VALUE.split("&")))
        counts = dict(summary.class_counts)
        self.assertEqual(counts[LINE], 8)   # L1M1 x5, L1MA2, L2a, X7D_LINE
        self.assertEqual(counts[SINE], 10)  # Alu* x6, MIRb x2, MIR, AmnSINE1
        self.assertEqual(counts[DNA], 5)    # Tigger1 x3, MADE1, Charlie26a
        # Sorted by count descending
        sorted_counts = [c for _, c in summary.class_counts]
        self.assertEqual(sorted_counts, sorted(sorted_counts, reverse=True))
