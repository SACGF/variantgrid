from django.conf import settings
from django.test import TestCase

from annotation.fake_annotation import get_fake_annotation_version
from library.genomics.vcf_enums import VCFSymbolicAllele
from snpdb.models import GenomeBuild, Variant, VariantCoordinate
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


class VariantTestCase(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.variant = slowly_create_test_variant("3", 128198980, 'A', 'T', cls.grch37)
        cls.ref_variant = slowly_create_test_variant("3", 128198980, 'A', 'A', cls.grch37)

        # Need this for HGVSMatcher in VariantCoordinate - that may be removed in future
        get_fake_annotation_version(cls.grch37)

    def test_variant_string(self):
        variant_string = "3:128198980 A>T"
        v = Variant.get_from_string(variant_string, self.grch37)
        self.assertEqual(self.variant, v)

    def test_variant_coordinate(self):
        vc = self.variant.coordinate
        v = Variant.get_from_variant_coordinate(vc, self.grch37)
        self.assertEqual(self.variant, v)

    def test_ref_variant_coordinate_explicit(self):
        vc_explicit = self.ref_variant.coordinate.as_external_explicit(self.grch37)
        self.assertEqual(vc_explicit.ref, vc_explicit.alt)

    def test_ref_variant_coordinate_internal(self):
        vc_symbolic = self.ref_variant.coordinate.as_internal_symbolic(self.grch37)
        self.assertEqual(vc_symbolic.alt, Variant.REFERENCE_ALT)

    def _test_internal_to_external_and_back(self, variant_coordinate, genome_build):
        internal = variant_coordinate.as_internal_symbolic(genome_build)
        external = variant_coordinate.as_external_explicit(genome_build)
        internal2 = external.as_internal_symbolic(genome_build)
        self.assertEqual(internal, internal2, msg="internal->external->internal")

        external2 = internal2.as_external_explicit(genome_build)
        self.assertEqual(external, external2, msg="external->internal->external")

    def _test_to_and_from_string(self, variant_coordinate, genome_build):
        variant_string = str(variant_coordinate)
        new_vc = VariantCoordinate.from_string(variant_string, genome_build)
        self.assertEqual(variant_coordinate, new_vc, msg="to/from string")

    def _test_coordinate_conversion(self, variant_coordinate, genome_build):
        self._test_internal_to_external_and_back(variant_coordinate, genome_build)
        self._test_to_and_from_string(variant_coordinate, genome_build)

    def test_del(self):
        vc = VariantCoordinate(chrom='21', position=47532744, ref='T', alt='<DEL>', svlen=-1224)
        self._test_coordinate_conversion(vc, self.grch37)

        vc2 = VariantCoordinate(chrom='7', position=23000505, ref='C', alt='<DEL>', svlen=-4084)
        self._test_coordinate_conversion(vc2, self.grch37)

    def test_dup(self):
        vc = VariantCoordinate(chrom='1', position=1000000, ref='T', alt='<DUP>', svlen=1000)
        self._test_coordinate_conversion(vc, self.grch37)

        vc2 = VariantCoordinate(chrom='7', position=41200841, ref='C', alt='<DUP>', svlen=2646)
        self._test_coordinate_conversion(vc2, self.grch37)

    def test_inv(self):
        vc = VariantCoordinate(chrom='10', position=89714001, ref='A', alt='<INV>', svlen=9549)
        self._test_coordinate_conversion(vc, self.grch37)

        vc = VariantCoordinate(chrom='1', position=100000, ref='C', alt='<INV>', svlen=1000)
        self._test_coordinate_conversion(vc, self.grch37)

    def test_multi_variant(self):
        ref = 'AACCGGTT'
        long_alt = "GATTACA" * 200  # 1.4kb
        vc = VariantCoordinate(chrom='10', position=89714001, ref=ref, alt=long_alt)
        vc_s = vc.as_internal_symbolic(self.grch37)
        # This is not symbolic, it's a multi-alt
        self.assertEquals(vc.ref, vc_s.ref)
        self.assertEquals(vc.alt, vc_s.alt)

    def test_delins(self):
        # clinvar variation 869248
        ref = 'AATATTTATATGCAGAGATATTGCTATTGCCTTAACCCAGAAATTATCACTGTTATTCTTTAGAATGGTGCAAAGAGGCATGATACATTGTATCATTATTGCCCTGAAAGAAAGAGATTAGGGAAAGTATTAGAAATAAGATAAACAAAAAAGTATATTAAAAGAAGAAAGCATTTTTTAAAATTACAAATGCAAAATTACCCTGATTTGGTCAATATGTGTACACATATTAAAACATTACACTTTAACCCATAAATATGTATAATGATTATGTATCAATTAAAAATAAAAGAAAATAAAGTAGGGAGATTATGAATATGCAAATAAGCACACATATATTCCAAATAGTAATGTACTAGGCAGACTGTGTAAAGTTTTTTTTTAAGTTACTTAATGTATCTCAGAGATATTTCCTTTTGTTATACACAATGTTAAGGCATTAAGTATAATAGTAAAAATTGCGGAGAAGAAAAAAAAAGAAAGCAAGAATTAAACAAAAGAAAACAATTGTTATGAACAGCAAATAAAAGAAACTAAAACGATCCTGAGACTTCCACACTGATGCAATCATTCGTCTGTTTCCCATTCTAAACTGTACCCTGTTACTTATCCCCTTCCTATGACATGAACTTAACCATAGAAAAGAAGGGGAAAGAAAACATCAAGCGTCCCATAGACTCACCCTGAAGTTCTCAGGATCCACGTGCAGCTTGTCACAGTGCAGCTCACTCAGTGTGGCAAAGGTGCCCTTGAGGTTGTCCAGGTGAGCCAGGCCATCACTAAAGGCACCGAGCACTTTCTTGCCATGAGCCTTCACCTTAGGGTTGCCCATAACAGCATCAGGAGTGGACAGATCCCCAAAGGACTCAAAGAACCTCTGGGTCCAAGGGTAGACCACCAGCAGCCTAAGGGTGGGAAAATAGACCAATAGGCAGAGAGAGTCAGTGCCTATCAGAAACCCAAGAGTCTTCTCTGTCTCCACATGCCCAGTTTCTATTGGTCTCCTTAAACCTGTCTTGTAACCTTGATACCAACCTGCCCAGGGCCTCACCACCAACTTCATCCACGTTCACCTTGCCCCACAGGGCAGTAACGGCAGACTTCTCCTCAGGAGTCAGATGCACCATGGTGTCTGTTTGAGGTTGCTAGTGAACACAGTTGTGTCAGAAGCAAATGTAAGCAATAGATGGCTCTGCCCTGACTTTTATGCCCAGCCCTGGCTCCTGCCCTCCCTGCTCCTGGGAGTAGATTGGCCAACCCTAGGGTGTGGCTCCACAGGGTGAGGTCTAAGTGATGACAGCCGTACCTGTCCTTGGCTCTTCTGGCACTGGCTTAGGAGTTGGACTTCAAACCCTCAGCCCTCCCTCTAAGATATATCTCTTGGCCCCATACCATCAGTACAAATTGCTACTAAAAACATCCTCCTTTGCAAGTGTATTTACGTAATATTTGGAATCACAGCTTGGTAAGCATATTGAAGATCGTTTTCCCAATTTTCTTATTACACAAATAAGAAG'
        vc = VariantCoordinate(chrom="11", position=5247125, ref=ref, alt="T")
        self._test_coordinate_conversion(vc, self.grch37)

    def test_one_exactly_symbolic_alt_size(self):
        alt = 'AGACAGAGTCTCGCTCTGTCGCCCAGGCTGGAGTGCAGTGCACAATCTTGGCTCACTGCAAGCTCCGCCTCCCAGGTTCACACCATTCTCCTGCCTCAGCCTCCCGAGTAGCCGGGACTACAGGCGCCCACCACCACGCCCAGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCATGTTAGTTAGCCAGGATGGTCTCGATCTCCTGACCTCGTGATCCACCCACCTCGGCCTCCCAAAGCACTGGGATTACAGGCATGAGCCACCGCGCCGAGCCCCAAGACCTTTCTTTATTACCAGGGCTTCCACAGACCTGACACATGGTAGTTCCTCAATAAATAATTGCAGAATTACTGAAAAATTTTACTGTTAACTTAGGCAGTGGTAAAACCATTGTTTGGTAGCTCAGAACTCAGCAAGTAAATAGCAACATTTGCTGGAAGAACAGATAGTTTTTCAAATCCAATTCAAGGACTGGGTATGGTGGCTCATGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCAGGCGTATCCAGGAGTTCGAGACTAGCCTGACCAACATGGTGAAACTCCGTCTCTACTAAAAATACAAAATTAGCCAGGTGTGGTGGTGGGCACCTGTAATCTCAGCTACTTGGGAGGCTGAGGCAGGAGAATCGCTTGAACCTGGTAGGCGGAGGTTGTAGTGAGCTGAGATTGTGCCATTGCTCTCCAGCCTGGGAAACAAGAGCAAAACTCCGTCTCAAAAAAAAAAAAAATCCAATTCAAATGATTATGGAAGTAGTGGAGAAATAAACAGGAAAATGATAAATAATTAAGATAATATATAATATGGCTATATTTTAATCTATTGTTGATATGATTTTCTCTTTTCCCCTTGGGATTAGTATCTATCTCTCTACTGGATATTAATTTGTTATATTTTCTCATTAGAGCAAGTTACTCAGATGGAAAACTGAAAGCCCCTCCTAAACCATGTGCTGGCAATCAAGG'
        vc = VariantCoordinate(chrom='3', position=37047542, ref='A', alt=alt)
        vc_symbolic = vc.as_internal_symbolic(self.grch37)
        self.assertEqual(len(alt), settings.VARIANT_SYMBOLIC_ALT_SIZE)  # If this fails, the test below won't work
        self.assertEqual(vc_symbolic.alt, VCFSymbolicAllele.DUP)

    def test_short_symbolic_alt(self):
        vc = VariantCoordinate(chrom='3', position=37047542, ref='A', alt="<DUP>",
                               svlen=settings.VARIANT_SYMBOLIC_ALT_SIZE-1)
        vc = vc.as_symbolic_or_explicit_according_to_size(self.grch37)
        self.assertIsNone(vc.svlen)
