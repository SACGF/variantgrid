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
        cls.chr21_gain = VariantCoordinate(chrom="21", position=1, ref="N", alt=VCFSymbolicAllele.CHROM_GAIN,
                                           svlen=GRCH38_CHR21_LENGTH - 1)

    def _canonical(self, chrom, position, ref, alt, svlen) -> VariantCoordinate:
        vc = VariantCoordinate(chrom=chrom, position=position, ref=ref, alt=alt, svlen=svlen)
        return vc.as_internal_canonical_form(self.grch38)

    def test_rule(self):
        chrx_length = 156040895
        for name, (chrom, position, ref, alt, svlen), expected_alt in [
            ("DRAGEN ploidy record", ("21", 1, "N", "<DUP>", GRCH38_CHR21_LENGTH - 1), VCFSymbolicAllele.CHROM_GAIN),
            ("From the first mappable base", ("21", 13_000_000, "A", "<DUP>", GRCH38_CHR21_LENGTH - 13_000_000),
             VCFSymbolicAllele.CHROM_GAIN),
            ("Loss", ("X", 1, "N", "<DEL>", -(chrx_length - 1)), VCFSymbolicAllele.CHROM_LOSS),
            ("Ends 18 Mb short of chr2", ("2", 1000, "A", "<DUP>", 223_000_000), VCFSymbolicAllele.DUP),
            ("No direction", ("21", 1, "N", "<CNV>", GRCH38_CHR21_LENGTH - 1), VCFSymbolicAllele.CNV),
            ("No cytobands", ("MT", 1, "G", "<DEL>", -16000), VCFSymbolicAllele.DEL),
        ]:
            with self.subTest(name):
                self.assertEqual(expected_alt, self._canonical(chrom, position, ref, alt, svlen).alt)

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
        vc = VariantCoordinate(chrom="21", position=1, ref="N", alt=VCFSymbolicAllele.CHROM_GAIN,
                               svlen=GRCH38_CHR21_LENGTH - 1)
        cls.allele = create_mock_allele(slowly_create_test_variant_from_coordinate(vc, GenomeBuild.grch38()),
                                        GenomeBuild.grch38())

    def test_same_chromosome(self):
        grch37 = GenomeBuild.grch37()
        expected = VariantCoordinate(chrom="21", position=1, ref="N", alt=VCFSymbolicAllele.CHROM_GAIN,
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
