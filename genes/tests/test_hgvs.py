"""
    PyHGVS has its own testing, this is specific to our code.
"""
from unittest import skip

from django.test import TestCase, override_settings
from pyhgvs import HGVSName  # This is used for pyhgvs specific test

from annotation.fake_annotation import get_fake_annotation_version
from annotation.tests.test_data_fake_genes import create_fake_transcript_version, create_gata2_transcript_version
from genes.hgvs import HGVSMatcher, HGVSException, HGVSConverterType
from genes.hgvs.hgvs_matcher import FakeTranscriptVersion
from genes.hgvs.pyhgvs.hgvs_converter_pyhgvs import PyHGVSVariant
from snpdb.models import GenomeBuild, VariantCoordinate


@override_settings(HGVS_VALIDATE_REFSEQ_TRANSCRIPT_LENGTH=False)
class TestHGVS(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(grch37)

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
            "NM_001754.5:c557T>A",  # Missing "."
            "NM_001754.5(RUNX1):557T>A",  # Has gene, missing "c."
            "NM_001754.5(RUNX1):c557T>A",  # Has gene, Missing "."
            "NC_000007.13:117199563G>T",  # Missing "g."
            "NC_000007.13:g117199563G>T",  # Missing "."
            "NM_000350.2(ABCA4):c-52delC",  # Missing "." with "-"
        ]

        hgvs_matcher = HGVSMatcher(genome_build=GenomeBuild.grch38(),
                                   hgvs_converter_type=HGVSConverterType.PYHGVS)
        for bad_hgvs in BAD_HGVS:
            try:
                hgvs_matcher.create_hgvs_variant(bad_hgvs)
                self.fail(f"Expected '{bad_hgvs}' to fail!")
            except:
                pass  # Expected

            fixed_hgvs = hgvs_matcher.clean_hgvs(bad_hgvs)[0]
            hgvs_matcher.create_hgvs_variant(fixed_hgvs)

    def test_fix_gene_transcript(self):
        swap_warning = "Swapped gene/transcript"
        uc_warning = "Upper cased transcript"

        TEST_CASES = [
            ("nm_000059.4:c.316+5G>A", [uc_warning]),
            ("nm_000059.4(BRCA1):c.316+5G>A", [uc_warning]),
            ("BRCA1(NM_000059.4):c.316+5G>A", [swap_warning]),
            ("BRCA1(nm_000059.4):c.316+5G>A",  [swap_warning, uc_warning]),
        ]
        hgvs_matcher = HGVSMatcher(GenomeBuild.grch38())
        for hgvs_string, expected_warnings in TEST_CASES:
            _, fix_messages = hgvs_matcher.fix_gene_transcript(hgvs_string)
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
            hgvs_variant = PyHGVSVariant(HGVSName(hgvs_string))
            hgvs_actual_trimmed = hgvs_variant.format(max_ref_length=10)
            self.assertEqual(hgvs_actual_trimmed, hgvs_expected_trimmed)

            hgvs_actual_no_trim = hgvs_variant.format(max_ref_length=100)
            self.assertEqual(hgvs_actual_no_trim, hgvs_string)  # No change

    @skip
    def test_c_hgvs_out_of_range(self):
        """ Disabled as it needs Ensembl contigs """
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        create_fake_transcript_version(genome_build)  # So we can lookup 'ENST00000300305.3'
        matcher = HGVSMatcher(genome_build)
        matcher.get_variant_tuple("ENST00000300305.3:c.1440A>T")

        with self.assertRaises(HGVSException):
            matcher.get_variant_tuple("ENST00000300305.3:c.9999A>T")

    def test_sort_transcript_versions(self):
        transcript_version_and_methods = [
            (FakeTranscriptVersion("", 1), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 1), HGVSMatcher.HGVS_METHOD_INTERNAL_LIBRARY),
            (FakeTranscriptVersion("", 2), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 2), HGVSMatcher.HGVS_METHOD_INTERNAL_LIBRARY),
            (FakeTranscriptVersion("", 3), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 3), HGVSMatcher.HGVS_METHOD_INTERNAL_LIBRARY),
            # Missing v4
            (FakeTranscriptVersion("", 4), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 5), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 5), HGVSMatcher.HGVS_METHOD_INTERNAL_LIBRARY),
            (FakeTranscriptVersion("", 6), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 6), HGVSMatcher.HGVS_METHOD_INTERNAL_LIBRARY),
        ]

        expected_up_then_down = [
            (FakeTranscriptVersion("", 5), HGVSMatcher.HGVS_METHOD_INTERNAL_LIBRARY),
            (FakeTranscriptVersion("", 6), HGVSMatcher.HGVS_METHOD_INTERNAL_LIBRARY),
            (FakeTranscriptVersion("", 3), HGVSMatcher.HGVS_METHOD_INTERNAL_LIBRARY),
            (FakeTranscriptVersion("", 2), HGVSMatcher.HGVS_METHOD_INTERNAL_LIBRARY),
            (FakeTranscriptVersion("", 1), HGVSMatcher.HGVS_METHOD_INTERNAL_LIBRARY),
            (FakeTranscriptVersion("", 4), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 5), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 6), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 3), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 2), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
            (FakeTranscriptVersion("", 1), HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY),
        ]

        expected_closest = [
            (FakeTranscriptVersion("", 5), HGVSMatcher.HGVS_METHOD_INTERNAL_LIBRARY),
            (FakeTranscriptVersion("", 3), HGVSMatcher.HGVS_METHOD_INTERNAL_LIBRARY),
            (FakeTranscriptVersion("", 6), HGVSMatcher.HGVS_METHOD_INTERNAL_LIBRARY),
            (FakeTranscriptVersion("", 2), HGVSMatcher.HGVS_METHOD_INTERNAL_LIBRARY),
            (FakeTranscriptVersion("", 1), HGVSMatcher.HGVS_METHOD_INTERNAL_LIBRARY),
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

    def test_hgvs_pyhgvs(self):
        self._test_hgvs_conversion(HGVSConverterType.PYHGVS)

    def test_hgvs_biocommons(self):
        # PyHGVS doesn't support INV but biocommons does
        extra_hgvs = [
            "NM_001145661.2(GATA2):c.1117_1131inv",
        ]
        self._test_hgvs_conversion(HGVSConverterType.BIOCOMMONS_HGVS, extra_hgvs)

    def _test_hgvs_conversion(self, hgvs_converter_type: HGVSConverterType, extra_hgvs=None):
        if extra_hgvs is None:
            extra_hgvs = []

        # GATA2 ClinVar
        HGVS_EXAMPLES = [
            "NM_001145661.2(GATA2):c.1121G>A",
            "NM_001145661.2(GATA2):c.1114G>A",
            "NM_001145661.2(GATA2):c.1017+572C>T",
            "NM_001145661.2(GATA2):c.1144-1G>C",
            # Single base del
            "NM_001145661.2(GATA2):c.1142del",
            "NM_001145661.2(GATA2):c.1124del",
            "NM_001145661.2(GATA2):c.1113del",
            # Multi-base del
            "NM_001145661.2(GATA2):c.1143+200_1198del",
            "NM_001145661.2(GATA2):c.1117_1131del",
            "NM_001145661.2(GATA2):c.1084_1095del",
            "NM_001145661.2(GATA2):c.1066_1095del",
            "NM_001145661.2(GATA2):c.1031_1049del",
            "NM_001145661.2(GATA2):c.1172_1175del",
            # "NM_001145661.2(GATA2):c.1017+513_1017+540del"
            # Biocommons converts this to 'NM_001145661.2(GATA2):c.1017+510_1017+537del'
            # clingen allele registry agrees with 'NM_001145661.2:c.1017+513_1017+540del'
            # This is the same coords, just not normalized as
            # HGVSUnsupportedOperationError: Normalization of intronic variants is not supported

            # Ins
            "NM_001145661.2(GATA2):c.1035_1036insTCTGGCC",
            "NM_001145661.2(GATA2):c.1034_1035insTCTTCTTGTGGCGGCTCTTCTGGCGGC",
            "NM_001145661.2(GATA2):c.1019_1020insCGACTGGGAGGGCAAGGCAG",
            # Delins
            "NM_001145661.2(GATA2):c.554_628delinsTAGCACCACGGGGGCT",
            "NM_001145661.2(GATA2):c.932_937delinsG",
            "NM_001145661.2(GATA2):c.405_409delinsGTA",
            "NM_001145661.2(GATA2):c.243delinsGC",
            # Dup
            "NM_001145661.2(GATA2):c.890_903dup",
            "NM_001145661.2(GATA2):c.1200_1216dup",
            "NM_001145661.2(GATA2):c.1126_1133dup",
            "NM_001145661.2(GATA2):c.1023_1038dup",
        ]

        genome_build = GenomeBuild.grch37()
        create_gata2_transcript_version(genome_build)
        matcher = HGVSMatcher(genome_build, hgvs_converter_type=hgvs_converter_type)
        for hgvs_string in HGVS_EXAMPLES + extra_hgvs:
            transcript_accession = matcher.get_transcript_accession(hgvs_string)
            vc = matcher.get_variant_tuple(hgvs_string)
            hgvs_variant = matcher.variant_coordinate_to_hgvs_variant(vc, transcript_accession)
            hgvs_out = hgvs_variant.format(max_ref_length=0)
            self.assertEqual(hgvs_string, hgvs_out, f"{hgvs_converter_type} Converting to and back to VariantCoordinate")

    def test_pyhgvs_reference_diff(self):
        self._test_reference_diff(HGVSConverterType.PYHGVS)

    def test_biocommons_reference_diff(self):
        self._test_reference_diff(HGVSConverterType.BIOCOMMONS_HGVS)

    def _test_reference_diff(self, hgvs_converter_type: HGVSConverterType):
        TEST_HGVS = [
            # No provided ref
            ('NM_001145661.2(GATA2):c.1113dup', VariantCoordinate(chrom='3', pos=128200691, ref='C', alt='CG'), True),
            # provided ref = genomic ref
            ('NM_001145661.2(GATA2):c.1113dupC', VariantCoordinate(chrom='3', pos=128200691, ref='C', alt='CG'), True),
            # GATA2 is '-' strand, so provided ref C = genomic ref = G
            # This is normalized left
            #('NM_001145661.2(GATA2):c.1113dupG', VariantCoordinate(chrom='3', pos=128200690, ref='G', alt='GC'),
            # HgvsMatchRefAllele(provided_ref='C', calculated_ref='G')),
        ]
        genome_build = GenomeBuild.grch37()
        create_gata2_transcript_version(genome_build)
        matcher = HGVSMatcher(genome_build, hgvs_converter_type=hgvs_converter_type)
        for hgvs_string, expected_vc, expected_matches_ref in TEST_HGVS:
            vcd = matcher.get_variant_tuple_used_transcript_kind_method_and_matches_reference(hgvs_string)
            #print(f"{vcd=}")
            #print(f"{type(vcd.matches_reference)=}")
            self.assertEqual(vcd.variant_coordinate, expected_vc)
            self.assertEqual(vcd.matches_reference, expected_matches_ref)
