from io import StringIO

from django.core.management import call_command
from django.test import TestCase

from library.genomics.vcf_enums import VCFSymbolicAllele
from snpdb.models import GenomeBuild, Variant, VariantCoordinate
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant_from_coordinate

GRCH38_CHR21_LENGTH = 46709983
GRCH38_CHR22_LENGTH = 50818468


class OneOffWholeChromosomeVariantsTest(TestCase):
    """ <DUP>s imported before the canonicaliser knew a whole chromosome (#2908) """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        grch38 = GenomeBuild.grch38()

        def stored_dup(chrom, position, svlen) -> Variant:
            vc = VariantCoordinate(chrom=chrom, position=position, ref="N", alt=VCFSymbolicAllele.DUP, svlen=svlen)
            return slowly_create_test_variant_from_coordinate(vc, grch38)

        cls.chr21_dup = stored_dup("21", 1, GRCH38_CHR21_LENGTH - 1)
        cls.focal_dup = stored_dup("21", 1, GRCH38_CHR21_LENGTH // 2 + 1)  # A candidate, not the whole chromosome
        cls.chr22_dup = stored_dup("22", 1, GRCH38_CHR22_LENGTH - 1)
        cls.chr22_gain = slowly_create_test_variant_from_coordinate(
            VariantCoordinate(chrom="22", position=1, ref="N", alt=VCFSymbolicAllele.DUP_CHROM,
                              svlen=GRCH38_CHR22_LENGTH - 1), grch38)

    def _run(self, **kwargs) -> str:
        out = StringIO()
        call_command("one_off_whole_chromosome_variants", stdout=out, **kwargs)
        return out.getvalue()

    def _coordinate(self, variant):
        return Variant.objects.get(pk=variant.pk).coordinate.as_tuple

    def test_dry_run_changes_nothing(self):
        before = self._coordinate(self.chr21_dup)
        output = self._run(dry_run=True)
        self.assertIn("1 variants to move", output)
        self.assertEqual(before, self._coordinate(self.chr21_dup))

    def test_moves_in_place(self):
        focal_before = self._coordinate(self.focal_dup)
        chr22_before = self._coordinate(self.chr22_dup)
        output = self._run()

        self.assertEqual(("21", 1, "N", VCFSymbolicAllele.DUP_CHROM, GRCH38_CHR21_LENGTH - 1),
                         self._coordinate(self.chr21_dup))
        self.assertEqual(GRCH38_CHR21_LENGTH, Variant.objects.get(pk=self.chr21_dup.pk).end)
        self.assertEqual(focal_before, self._coordinate(self.focal_dup))
        self.assertEqual(chr22_before, self._coordinate(self.chr22_dup), "Its whole-chromosome form already exists")
        self.assertIn(f"already exists as variant {self.chr22_gain.pk}", output)
