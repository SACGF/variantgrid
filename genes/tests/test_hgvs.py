"""
    PyHGVS has its own testing, this is specific to our code.
"""
from unittest import skip

from django.test.testcases import TestCase
from pyhgvs import HGVSName, InvalidHGVSName

from annotation.tests.test_data_fake_genes import create_fake_transcript_version
from genes.hgvs import HGVSMatcher, FakeTranscriptVersion, HGVSNameExtra
from snpdb.models import GenomeBuild


class TestHGVS(TestCase):

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
            # Unbalanced brackets
            "NM_001754.5):c.557T>A",
            "(NM_004991.4:c.2577+4A>T",
            # Good brackets HGVS (just testing gene symbol)
            "NM_001754.5(RUNX1):c.1415T>C",
            "NM_032638:c.1126_1133DUP",  # Case
            "NM_001754.5:557T>A",  # Missing "c."
            "NC_000007.13:117199563G>T",  # Missing "g."
        ]

        for bad_hgvs in BAD_HGVS:
            try:
                HGVSName(bad_hgvs)
                self.fail(f"Expected '{bad_hgvs}' to fail!")
            except:
                pass  # Expected

            fixed_hgvs = HGVSMatcher.clean_hgvs(bad_hgvs)[0]
            HGVSName(fixed_hgvs)

    def test_fix_gene_transcript(self):
        swap_warning = "Warning: swapped gene/transcript"
        uc_warning = "Warning: Upper cased transcript"

        TEST_CASES = [
            ("nm_000059.4:c.316+5G>A", [uc_warning]),
            ("nm_000059.4(BRCA1):c.316+5G>A", [uc_warning]),
            ("BRCA1(NM_000059.4):c.316+5G>A", [swap_warning]),
            ("BRCA1(nm_000059.4):c.316+5G>A",  [swap_warning, uc_warning]),
        ]
        for hgvs_string, expected_warnings in TEST_CASES:
            _, fix_messages = HGVSMatcher.fix_gene_transcript(hgvs_string)
            for ew in expected_warnings:
                self.assertIn(ew, fix_messages, f"Warning for {hgvs_string}")

    def test_format_hgvs_remove_long_ref(self):
        LONG_AND_TRIMMED_HGVS = {  # 10bp
            "NM_000726.4(CACNB4):c.162_173delCTACACAAGCAGinsGA": "NM_000726.4(CACNB4):c.162_173delinsGA",
            # "NM_007294.3:c.5080_5090del11insAA": "NM_007294.3:c.5080_5090delinsAA",
            "NM_007294.3:c.5080_5090delCCCCCCCCCCCinsAA": "NM_007294.3:c.5080_5090delinsAA",
            "NM_000726.4(CACNB4):c.162_173delCTACACAAGCAG": "NM_000726.4(CACNB4):c.162_173del",
        }

        for hgvs_string, hgvs_expected_trimmed in LONG_AND_TRIMMED_HGVS.items():
            hgvs_name_extra = HGVSNameExtra(HGVSName(hgvs_string))
            hgvs_actual_trimmed = hgvs_name_extra.format(max_ref_length=10)
            self.assertEqual(hgvs_actual_trimmed, hgvs_expected_trimmed)

            hgvs_actual_no_trim = hgvs_name_extra.format(max_ref_length=100)
            self.assertEqual(hgvs_actual_no_trim, hgvs_string)  # No change

    @skip
    def test_c_hgvs_out_of_range(self):
        """ Disabled as it needs Ensembl contigs """
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        create_fake_transcript_version(genome_build)  # So we can lookup 'ENST00000300305.3'
        matcher = HGVSMatcher(genome_build)
        matcher.get_variant_tuple("ENST00000300305.3:c.1440A>T")

        with self.assertRaises(InvalidHGVSName):
            matcher.get_variant_tuple("ENST00000300305.3:c.9999A>T")

    def test_sort_transcript_versions(self):
        transcript_version_and_methods = [
            (FakeTranscriptVersion("", 1), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 1), HGVSMatcher.HGVS_METHOD_PYHGVS),
            (FakeTranscriptVersion("", 2), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 2), HGVSMatcher.HGVS_METHOD_PYHGVS),
            (FakeTranscriptVersion("", 3), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 3), HGVSMatcher.HGVS_METHOD_PYHGVS),
            # Missing v4
            (FakeTranscriptVersion("", 4), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 5), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 5), HGVSMatcher.HGVS_METHOD_PYHGVS),
            (FakeTranscriptVersion("", 6), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 6), HGVSMatcher.HGVS_METHOD_PYHGVS),
        ]

        expected_up_then_down = [
            (FakeTranscriptVersion("", 5), HGVSMatcher.HGVS_METHOD_PYHGVS),
            (FakeTranscriptVersion("", 6), HGVSMatcher.HGVS_METHOD_PYHGVS),
            (FakeTranscriptVersion("", 3), HGVSMatcher.HGVS_METHOD_PYHGVS),
            (FakeTranscriptVersion("", 2), HGVSMatcher.HGVS_METHOD_PYHGVS),
            (FakeTranscriptVersion("", 1), HGVSMatcher.HGVS_METHOD_PYHGVS),
            (FakeTranscriptVersion("", 4), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 5), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 6), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 3), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 2), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 1), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
        ]

        expected_closest = [
            (FakeTranscriptVersion("", 5), HGVSMatcher.HGVS_METHOD_PYHGVS),
            (FakeTranscriptVersion("", 3), HGVSMatcher.HGVS_METHOD_PYHGVS),
            (FakeTranscriptVersion("", 6), HGVSMatcher.HGVS_METHOD_PYHGVS),
            (FakeTranscriptVersion("", 2), HGVSMatcher.HGVS_METHOD_PYHGVS),
            (FakeTranscriptVersion("", 1), HGVSMatcher.HGVS_METHOD_PYHGVS),
            (FakeTranscriptVersion("", 4), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 5), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 3), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 6), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 2), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 1), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
        ]

        version = 4
        key_up_then_down = HGVSMatcher._get_sort_key_transcript_version_and_methods(version)
        sorted_up_then_down = list(sorted(transcript_version_and_methods, key=key_up_then_down))
        self.assertEqual(sorted_up_then_down, expected_up_then_down, "Sorted up then down")

        key_closest = HGVSMatcher._get_sort_key_transcript_version_and_methods(version, closest=True)
        sorted_closest = list(sorted(transcript_version_and_methods, key=key_closest))
        self.assertEqual(sorted_closest, expected_closest, "Sorted closest")
