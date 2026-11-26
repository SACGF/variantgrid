from unittest.mock import patch

from django.test import TestCase

from library.genomics.vcf_enums import VCFSymbolicAllele
from snpdb.enums import AlleleConversionTool
from snpdb.liftover import _liftover_using_dest_variant_coordinate
from snpdb.models import GenomeBuild, Variant, VariantCoordinate
from snpdb.tests.utils.vcf_testing_utils import (
    create_mock_allele,
    slowly_create_test_variant_from_coordinate,
)

GRCH38_CHR21_LENGTH = 46709983
GRCH37_CHR21_LENGTH = 48129895


class WholeChromosomeCanonicalFormTest(TestCase):
    """ #2908 - against the cytobands the migrations load """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch38 = GenomeBuild.grch38()
        cls.chr21_gain = VariantCoordinate(chrom="21", position=1, ref="N", alt=VCFSymbolicAllele.DUP_CHROM,
                                           svlen=GRCH38_CHR21_LENGTH - 1)

    def _canonical(self, chrom, position, ref, alt, svlen) -> VariantCoordinate:
        vc = VariantCoordinate(chrom=chrom, position=position, ref=ref, alt=alt, svlen=svlen)
        return vc.as_internal_canonical_form(self.grch38)

    def test_rule(self):
        """ The span callers cover, to within a 1 Mb slack each end - and the near misses an import reports """
        chrx_length = 156040895
        dup, dup_chrom = VCFSymbolicAllele.DUP, VCFSymbolicAllele.DUP_CHROM
        for name, (chrom, position, alt, end), expected_alt, near_miss in [
            ("DRAGEN ploidy record", ("21", 1, "<DUP>", GRCH38_CHR21_LENGTH), dup_chrom, False),
            ("From chr21's first non-N base", ("21", 5_010_001, "<DUP>", GRCH38_CHR21_LENGTH), dup_chrom, False),
            ("From the start of 21q", ("21", 13_000_000, "<DUP>", GRCH38_CHR21_LENGTH - 10_000), dup_chrom, False),
            ("Stops short of 21q22.3", ("21", 13_000_000, "<DUP>", 41_100_000), dup, False),
            ("Long arm", ("1", 143_200_001, "<DUP>", 248_946_422), dup, False),
            ("chrY between its PARs", ("Y", 2_781_480, "<DUP>", 26_600_000), dup_chrom, False),
            ("Loss", ("X", 1, "<DEL>", chrx_length), VCFSymbolicAllele.DEL_CHROM, False),
            ("Ends 18 Mb short of chr2", ("2", 1000, "<DUP>", 223_001_000), dup, True),
            ("No direction", ("21", 1, "<CNV>", GRCH38_CHR21_LENGTH), VCFSymbolicAllele.CNV, False),
            ("No cytobands", ("MT", 1, "<DEL>", 16001), VCFSymbolicAllele.DEL, False),
        ]:
            with self.subTest(name):
                svlen = end - position
                if alt == "<DEL>":
                    svlen = -svlen
                vc = VariantCoordinate(chrom=chrom, position=position, ref="N", alt=alt, svlen=svlen)
                self.assertEqual(expected_alt, vc.as_internal_canonical_form(self.grch38).alt)
                self.assertEqual(near_miss, vc.is_whole_chromosome_near_miss(self.grch38))

    def test_one_representation(self):
        from_13mb = self._canonical("21", 13_000_000, "A", "<DUP>", GRCH38_CHR21_LENGTH - 13_000_000)
        self.assertEqual(self.chr21_gain, from_13mb)
        self.assertEqual(GRCH38_CHR21_LENGTH, from_13mb.end)
        self.assertEqual(self.chr21_gain, self.chr21_gain.as_internal_canonical_form(self.grch38))
        self.assertEqual(self.chr21_gain, VariantCoordinate.from_string(str(self.chr21_gain), self.grch38))

    def test_no_hgvs_or_clingen(self):
        variant = slowly_create_test_variant_from_coordinate(self.chr21_gain, self.grch38)
        self.assertTrue(variant.is_whole_chromosome)
        self.assertEqual([variant], list(Variant.objects.filter(Variant.get_whole_chromosome_q())))
        self.assertIsNone(self.chr21_gain.symbolic_hgvs_interval)
        self.assertFalse(variant.can_make_g_hgvs)
        self.assertFalse(variant.can_have_c_hgvs)
        self.assertFalse(variant.can_have_clingen_allele)


class WholeChromosomeLiftoverTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        vc = VariantCoordinate(chrom="21", position=1, ref="N", alt=VCFSymbolicAllele.DUP_CHROM,
                               svlen=GRCH38_CHR21_LENGTH - 1)
        cls.allele = create_mock_allele(slowly_create_test_variant_from_coordinate(vc, GenomeBuild.grch38()),
                                        GenomeBuild.grch38())

    def test_same_chromosome(self):
        grch37 = GenomeBuild.grch37()
        expected = VariantCoordinate(chrom="21", position=1, ref="N", alt=VCFSymbolicAllele.DUP_CHROM,
                                     svlen=GRCH37_CHR21_LENGTH - 1)
        self.assertEqual([(AlleleConversionTool.WHOLE_CHROMOSOME, expected, None)],
                         list(_liftover_using_dest_variant_coordinate(self.allele, grch37, failed_tools=set())))

    def test_destination_without_the_chromosome(self):
        grch37 = GenomeBuild.objects.get(name="GRCh37")  # Its own chrom_contig_mappings
        with patch.object(GenomeBuild, "chrom_contig_mappings", {}):
            [(conversion_tool, variant_coordinate, error)] = _liftover_using_dest_variant_coordinate(
                self.allele, grch37, failed_tools=set())
        self.assertEqual(AlleleConversionTool.WHOLE_CHROMOSOME, conversion_tool)
        self.assertIsNone(variant_coordinate)
        self.assertIn("no chromosome '21'", error)
