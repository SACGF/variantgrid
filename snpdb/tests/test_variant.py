from bioutils.sequences import reverse_complement
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
        cls.deletion_variant = slowly_create_test_variant("3", 128198990, 'GAG', 'G', cls.grch37)
        cls.insertion_variant = slowly_create_test_variant("3", 128198995, 'G', 'GAG', cls.grch37)

        cls.contig_1 = cls.grch37.contigs.get(name="1")

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

    def test_internal_representation_ref(self):
        vc_ref = VariantCoordinate(chrom='21', position=47532744, ref='T', alt='T')
        vc_internal = vc_ref.as_internal_canonical_form(self.grch37)
        self.assertEqual(vc_internal.alt, Variant.REFERENCE_ALT)

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
        self.assertEqual(vc.ref, vc_s.ref)
        self.assertEqual(vc.alt, vc_s.alt)

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
        vc = vc.as_internal_canonical_form(self.grch37)
        self.assertIsNone(vc.svlen)

    def test_inv_at_exactly_symbolic_alt_size_stays_explicit(self):
        # INV intentionally uses strictly > (not >=) so an inversion of exactly
        # VARIANT_SYMBOLIC_ALT_SIZE stays as an explicit sequence.  Changing
        # the threshold would require re-checking/converting existing INV records.
        # "AC" repeats have reverse complement "GT" repeats — non-palindromic.
        size = settings.VARIANT_SYMBOLIC_ALT_SIZE
        ref = ("AC" * (size // 2 + 1))[:size]
        vc = VariantCoordinate(chrom="1", position=1_000_000, ref=ref, alt=reverse_complement(ref))
        vc_s = vc.as_internal_symbolic(self.grch37)
        self.assertNotEqual(vc_s.alt, VCFSymbolicAllele.INV,
            "INV threshold was deliberately left as strictly > (not >=). If you changed it to >=, "
            "existing production databases may contain explicit INV sequences of exactly "
            f"VARIANT_SYMBOLIC_ALT_SIZE={size} bp that will no longer round-trip correctly. "
            "You must write a migration to convert those records to <INV> symbolic form before "
            "changing this threshold.")

    def test_inv_above_symbolic_alt_size_becomes_symbolic(self):
        size = settings.VARIANT_SYMBOLIC_ALT_SIZE + 1
        ref = ("AC" * (size // 2 + 1))[:size]
        vc = VariantCoordinate(chrom="1", position=1_000_000, ref=ref, alt=reverse_complement(ref))
        vc_s = vc.as_internal_symbolic(self.grch37)
        self.assertEqual(vc_s.alt, VCFSymbolicAllele.INV)

    def test_as_internal_canonical_form(self):
        # test for issue #1214 - however we get there, it should end up the same
        vc_symbolic = VariantCoordinate(chrom='1', position=1000000, ref='T', alt='<DUP>', svlen=999)
        vc_explicit = vc_symbolic.as_external_explicit(self.grch37)

        vc_from_symbolic = vc_symbolic.as_internal_canonical_form(self.grch37)
        vc_from_explicit = vc_explicit.as_internal_canonical_form(self.grch37)
        self.assertEqual(vc_from_symbolic, vc_from_explicit)

    def test_is_deletion(self):
        """ Issue #1449 BUG 1 - len(self.locus.ref) was missing .seq, causing TypeError """
        self.assertTrue(self.deletion_variant.is_deletion)
        self.assertFalse(self.insertion_variant.is_deletion)
        self.assertFalse(self.variant.is_deletion)

    def test_is_insertion(self):
        self.assertTrue(self.insertion_variant.is_insertion)
        self.assertFalse(self.deletion_variant.is_insertion)
        self.assertFalse(self.variant.is_insertion)

    def test_clingen_allele_size(self):
        """ Issue #1449 BUG 2 - len(self.locus.ref) + len(self.alt) were missing .seq, causing TypeError """
        self.assertEqual(self.deletion_variant._clingen_allele_size, len('GAG') + len('G'))
        self.assertEqual(self.insertion_variant._clingen_allele_size, len('G') + len('GAG'))
        self.assertEqual(self.variant._clingen_allele_size, len('A') + len('T'))

    def test_can_have_c_hgvs(self):
        """ Issue #1449 BUG 3 - operator precedence caused TypeError for reference variants (svlen=None, can_have_annotation=False) """
        self.assertTrue(self.variant.can_have_c_hgvs)
        self.assertFalse(self.ref_variant.can_have_c_hgvs)

    def test_symbolic_from_string(self):
        variant_string_lower = "1:10000-50000 <del>"  # Written with lower case
        vc_lower = VariantCoordinate.from_string(variant_string_lower, self.grch37)
        vc_upper = VariantCoordinate.from_string(variant_string_lower.upper(), self.grch37)
        self.assertEqual(vc_lower, vc_upper, msg="symbolic from string case insensitive")

    # Variant.validate() — position boundary checks

    def test_validate_last_base_of_contig_is_valid(self):
        # Regression: validate() used `position < contig.length` (strict), rejecting the last base.
        last_base = self.contig_1.length
        errors = Variant.validate(self.grch37, "1", last_base)
        self.assertEqual(errors, [],
            f"Position {last_base} is valid but validate() returned: {errors}")

    def test_validate_position_beyond_contig_is_invalid(self):
        errors = Variant.validate(self.grch37, "1", self.contig_1.length + 1)
        self.assertGreater(len(errors), 0)

    # Variant Q-object classifiers vs instance properties

    def _in_qs(self, q, variant):
        return Variant.objects.filter(q, pk=variant.pk).exists()

    def test_snp_q_matches_snp(self):
        self.assertTrue(self._in_qs(Variant.get_snp_q(), self.variant))

    def test_snp_q_does_not_match_reference(self):
        self.assertFalse(self._in_qs(Variant.get_snp_q(), self.ref_variant))

    def test_is_deletion_property_agrees_with_q(self):
        for v in (self.deletion_variant, self.variant, self.insertion_variant, self.ref_variant):
            self.assertEqual(v.is_deletion, self._in_qs(Variant.get_deletion_q(), v),
                f"{v}: is_deletion disagrees with get_deletion_q()")

    def test_is_insertion_property_agrees_with_q(self):
        for v in (self.insertion_variant, self.variant, self.deletion_variant, self.ref_variant):
            self.assertEqual(v.is_insertion, self._in_qs(Variant.get_insertion_q(), v),
                f"{v}: is_insertion disagrees with get_insertion_q()")

    def test_is_reference_property_agrees_with_q(self):
        for v in (self.ref_variant, self.variant, self.deletion_variant, self.insertion_variant):
            self.assertEqual(v.is_reference, self._in_qs(Variant.get_reference_q(), v),
                f"{v}: is_reference disagrees with get_reference_q()")

    def test_indel_q_matches_deletion_and_insertion_not_snp(self):
        self.assertTrue(self._in_qs(Variant.get_indel_q(), self.deletion_variant))
        self.assertTrue(self._in_qs(Variant.get_indel_q(), self.insertion_variant))
        self.assertFalse(self._in_qs(Variant.get_indel_q(), self.variant))
