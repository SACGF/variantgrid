from django.test import TestCase, override_settings
from genes.hgvs import HGVSMatcher, HGVSConverterType
from snpdb.models import GenomeBuild

@override_settings(HGVS_VALIDATE_REFSEQ_TRANSCRIPT_LENGTH=False)
class TestCleanHGVS(TestCase):
    def setUp(self):
        self.hgvs_test_instance = HGVSMatcher(genome_build=GenomeBuild.grch38(),
                                   hgvs_converter_type=HGVSConverterType.PYHGVS)

    def test_clean_hgvs(self):
        test_cases = [
            ("M_206933.3(USH2A):c.4298G>A", "NM_206933.3(USH2A):c.4298G>A"),
            ("NM_206933.3(USH2A)::c.4298G>A", "NM_206933.3(USH2A):c.4298G>A"),
            ("M_206933.3USH2A):c.4298G>A", "NM_206933.3USH2A:c.4298G>A"),
            ("m_206933.3USH2A):c.4298G>A", "NM_206933.3USH2A:c.4298G>A"),
            ("nm_206933.3(USH2A):c.4298G>A", "NM_206933.3(USH2A):c.4298G>A"),
            ("m_206933.3(USH2A)::c.4298G>A", "NM_206933.3(USH2A):c.4298G>A"),
            ("m_206933.3(USH2A):c.4298G>A", "NM_206933.3(USH2A):c.4298G>A"),
        ]

        for input_hgvs, expected_result in test_cases:
            result = self.hgvs_test_instance.clean_hgvs(input_hgvs)
            self.assertEqual(result[0], expected_result)
