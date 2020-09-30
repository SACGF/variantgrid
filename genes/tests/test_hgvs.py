"""
    PyHGVS has its own testing, this is specific to our code.
"""
from django.test.testcases import TestCase
from pyhgvs import HGVSName
from genes.hgvs import HGVSMatcher


class TestAnnotationVCF(TestCase):

    def test_clean_hgvs(self):
        BAD_HGVS = [
            "NM_205768 c.44A>G",  # Missing colon (no version)
            "NM_005629.3:c1403A>C",  # Missing dot after kind
            "NM_001101.4 c.95C>G",  # Missing colon
            "NM_00380.3: c.648_649delGA",  # space after colon
            "NC_000023.10:g. 31496384G>A",
            "NM_004245: :c.337G>T",  # Double colon
            "NC_000017.10:g.21085664 G>C",  # Space after numbers
            "NC_000023.10:g. 133547943G>A",  # Space after g.
            # Missing transcript underscore, Missing colon, Missing dot after g
            # Space between position and reference base
            "NC000002.10g39139341 C>T",
        ]

        for bad_hgvs in BAD_HGVS:
            try:
                HGVSName(bad_hgvs)
                self.fail(f"Expected '{bad_hgvs}' to fail!")
            except:
                pass  # Expected

            fixed_hgvs = HGVSMatcher.clean_hgvs(bad_hgvs)
            HGVSName(fixed_hgvs)

    def test_format_hgvs_remove_long_ref(self):
        LONG_HGVS_OLD_NEW = {
            # "NM_000726.4(CACNB4):c.162_173delCTACACAAGCAGinsGA": "NM_000726.4(CACNB4):c.162_173delinsGA",
            # "NM_007294.3:c.5080_5090del11insAA": "NM_007294.3:c.5080_5090delinsAA",
            "NM_000726.4(CACNB4):c.162_173delCTACACAAGCAG": "NM_000726.4(CACNB4):c.162_173del",
        }

        for hgvs_string, hgvs_expected in LONG_HGVS_OLD_NEW.items():
            hgvs_name = HGVSName(hgvs_string)
            HGVSMatcher.format_hgvs_remove_long_ref(hgvs_name)
            self.assertEquals(hgvs_name.format(), hgvs_expected)
